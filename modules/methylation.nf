
process methylDackel_mbias {
    label 'medium_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts/mbias"
    // Add specific error strategy for this process to handle memory issues
    errorStrategy = { task.exitStatus in [143,137,104,134,139] || task.attempt <= 3 ? 'retry' : 'finish' }
    // Increase max retries for this process
    maxRetries = 3
    // Add memory directive to ensure adequate memory allocation
    memory = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def fileSizeGB = md_bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
        // Calculate memory based on file size with a higher multiplier
        def fileBasedMemGB = Math.ceil(fileSizeGB * 1.5).doubleValue() // Increased multiplier from 1.2 to 1.5
        def minMemGB = slurm_profile ? 4.GB : 2.GB // Increased minimum memory

        // Use the maximum of calculated memory or minimum memory, multiplied by attempt number
        def memToUse = Math.max(minMemGB.toGiga(), fileBasedMemGB).GB * task.attempt

        // Cap at max_memory
        check_max(memToUse, 'memory')
    }
    // Add time directive to ensure adequate time allocation
    time = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def default_time = slurm_profile ? 24.h : 8.h // Increased time allocation

        check_max(default_time * task.attempt, 'time')
    }
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"  // Skip procps-ng for macOS
        } else {
            "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8 conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }
    publishDir "${params.outputDir}/methylDackelExtracts/mbias"

    input:
        tuple val(library), path(md_bam), path(md_bai), val(barcodes)
        path(genome_path)

    output:
        path('*.svg'), emit: mbias_output_svg
        path('*.tsv'), emit: mbias_output_tsv
        tuple val(params.email), val(library), path('*.tsv'), emit: for_agg

    script:
    // Set memory and CPU resources based on BAM file size
    def fileSizeGB = md_bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
    def currentMemoryGB = task.memory.toGiga() // Convert task.memory to GB

    // Ensure all values are explicitly converted to double to avoid type ambiguity
    def maxMemoryGB = params.max_memory.toGiga().doubleValue()
    def currentMemGB = currentMemoryGB.doubleValue()
    def minMemGB = 1.0d
    def fileBasedMemGB = Math.ceil(fileSizeGB * 1.5).doubleValue() // Increased multiplier from 1.2 to 1.5

    // MethylDackel mbias is less memory-intensive than extract
    // but still benefits from scaling with file size
    def memoryGB = Math.min(
        maxMemoryGB,
        Math.max(Math.max(currentMemGB, minMemGB), fileBasedMemGB)
    )

    // Adjust CPUs based on available resources
    // Convert to integers for CPU calculations
    def taskCpusInt = task.cpus.intValue()
    def fileSizeInt = Math.ceil(fileSizeGB).intValue()
    def minCpus = 2
    def maxCpus = 4

    def cpusToUse = Math.min(
        taskCpusInt,
        Math.max(minCpus, Math.min(maxCpus, fileSizeInt))
    )

    task.memory = "${memoryGB} GB"

    """
    echo "Input BAM size: ${fileSizeGB} GB"
    echo "Memory allocated for this task: ${task.memory}"
    echo "CPUs allocated: ${cpusToUse}"

    # Monitor memory usage
    echo "Available memory before processing:"
    free -h || echo "free command not available"

    # Create a dedicated temp directory with random name to avoid conflicts
    TEMP_DIR="\${TMPDIR:-${params.tmp_dir}}/mbias_\${RANDOM}"
    mkdir -p "\$TEMP_DIR"
    echo "Using temporary directory: \$TEMP_DIR"

    # Check disk space in temp directory
    echo "Disk space in temp directory:"
    df -h "\$TEMP_DIR" || echo "df command not available"

    # Find the genome file
    genome=\$(ls *fa)

    # Step 1: Create the combined mbias TSV file
    echo "Step 1: Creating combined mbias TSV file..."
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ${library}_${barcodes}_combined_mbias.tsv

    # Use parallelization for samtools view
    echo "Getting chromosome list..."
    chrs=(`samtools view -@ ${cpusToUse} -H ${md_bam} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

    # Check if we got any chromosomes
    if [ \${#chrs[@]} -eq 0 ]; then
        echo "Warning: No chromosomes found in BAM header. Using first 10 reads to determine chromosome."
        chrs=(`samtools view -@ ${cpusToUse} ${md_bam} | head -10 | cut -f 3 | sort | uniq`)
        if [ \${#chrs[@]} -eq 0 ]; then
            echo "Error: Could not determine any chromosomes from BAM file."
            exit 1
        fi
    fi

    echo "Found \${#chrs[@]} chromosomes. Using first chromosome: \${chrs[0]}"

    # Process each chromosome and context
    for chr in \${chrs[*]}; do
        echo "Processing chromosome: \$chr"
        for context in CHH CHG CpG; do
            echo "Processing context: \$context"
            arg=''
            if [ "\$context" = 'CHH' ]; then
                arg='--CHH --noCpG'
            elif [ "\$context" = 'CHG' ]; then
                arg='--CHG --noCpG'
            fi

            # Run first MethylDackel mbias command
            echo "Running MethylDackel mbias for \$context without duplicates..."
            MethylDackel mbias --noSVG \$arg -@ ${cpusToUse} -r \$chr \${genome} ${md_bam} > "\$TEMP_DIR/\${chr}_\${context}_nodups.txt"

            mbias1_exit=\$?
            if [ \$mbias1_exit -ne 0 ]; then
                echo "MethylDackel mbias (no duplicates) failed with exit code \$mbias1_exit"
                echo "This might be due to memory constraints."

                # Check memory usage
                echo "Current memory usage:"
                free -h || echo "free command not available"

                # If this is not the last retry, exit with a code that will trigger a retry
                if [ ${task.attempt} -lt 3 ]; then
                    echo "Will retry with more memory"
                    exit 137  # Memory error code that will trigger retry
                fi

                exit \$mbias1_exit
            fi

            # Run second MethylDackel mbias command with duplicates
            echo "Running MethylDackel mbias for \$context with duplicates..."
            MethylDackel mbias --noSVG --keepDupes -F 2816 \$arg -@ ${cpusToUse} -r \$chr \${genome} ${md_bam} > "\$TEMP_DIR/\${chr}_\${context}_withdups.txt"

            mbias2_exit=\$?
            if [ \$mbias2_exit -ne 0 ]; then
                echo "MethylDackel mbias (with duplicates) failed with exit code \$mbias2_exit"
                echo "This might be due to memory constraints."

                # Check memory usage
                echo "Current memory usage:"
                free -h || echo "free command not available"

                # If this is not the last retry, exit with a code that will trigger a retry
                if [ ${task.attempt} -lt 3 ]; then
                    echo "Will retry with more memory"
                    exit 137  # Memory error code that will trigger retry
                fi

                exit \$mbias2_exit
            fi

            # Join the results
            echo "Joining results for \$chr \$context..."
            join -t \$'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                cat "\$TEMP_DIR/\${chr}_\${context}_nodups.txt" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            <( \
                cat "\$TEMP_DIR/\${chr}_\${context}_withdups.txt" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            | sed "s/^/\${chr}\t\${context}\t/" \
            >> ${library}_${barcodes}_combined_mbias.tsv
        done
    done

    # Verify the combined mbias file exists and has content
    if [ ! -s "${library}_${barcodes}_combined_mbias.tsv" ]; then
        echo "Error: Combined mbias TSV file is empty or does not exist"
        exit 1
    fi

    # Step 2: Generate SVG files for CHN contexts
    echo "Step 2: Generating SVG files for CHN contexts..."
    MethylDackel mbias -@ ${cpusToUse} --noCpG --CHH --CHG -r \${chrs[0]} \${genome} ${md_bam} ${library}_chn

    mbias_chn_exit=\$?
    if [ \$mbias_chn_exit -ne 0 ]; then
        echo "MethylDackel mbias for CHN contexts failed with exit code \$mbias_chn_exit"
        echo "This might be due to memory constraints."

        # Check memory usage
        echo "Current memory usage:"
        free -h || echo "free command not available"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry with more memory"
            exit 137  # Memory error code that will trigger retry
        fi

        exit \$mbias_chn_exit
    fi

    # Check OS type and use appropriate sed syntax for CHN files
    echo "Updating CHN SVG files..."
    if [[ "\$(uname)" == "Darwin" ]]; then
        # macOS
        for f in *chn*.svg; do 
            if [ -f "\$f" ]; then
                sed -i '' "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f
            else
                echo "Warning: Expected CHN SVG file \$f not found"
            fi
        done
    else
        # Linux and other Unix-like systems
        for f in *chn*.svg; do 
            if [ -f "\$f" ]; then
                sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f
            else
                echo "Warning: Expected CHN SVG file \$f not found"
            fi
        done
    fi

    # Step 3: Generate SVG files for CpG context
    echo "Step 3: Generating SVG files for CpG context..."
    MethylDackel mbias -@ ${cpusToUse} -r \${chrs[0]} \${genome} ${md_bam} ${library}_cpg

    mbias_cpg_exit=\$?
    if [ \$mbias_cpg_exit -ne 0 ]; then
        echo "MethylDackel mbias for CpG context failed with exit code \$mbias_cpg_exit"
        echo "This might be due to memory constraints."

        # Check memory usage
        echo "Current memory usage:"
        free -h || echo "free command not available"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry with more memory"
            exit 137  # Memory error code that will trigger retry
        fi

        exit \$mbias_cpg_exit
    fi

    # Check OS type and use appropriate sed syntax for CpG files
    echo "Updating CpG SVG files..."
    if [[ "\$(uname)" == "Darwin" ]]; then
        # macOS
        for f in *cpg*.svg; do 
            if [ -f "\$f" ]; then
                sed -i '' "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f
            else
                echo "Warning: Expected CpG SVG file \$f not found"
            fi
        done
    else
        # Linux and other Unix-like systems
        for f in *cpg*.svg; do 
            if [ -f "\$f" ]; then
                sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f
            else
                echo "Warning: Expected CpG SVG file \$f not found"
            fi
        done
    fi

    # Verify SVG files exist
    svg_count=\$(ls -1 *.svg 2>/dev/null | wc -l)
    if [ \$svg_count -eq 0 ]; then
        echo "Warning: No SVG files were generated"
    else
        echo "Generated \$svg_count SVG files"
    fi

    # Clean up temp directory
    rm -rf "\$TEMP_DIR"

    echo "MethylDackel mbias process completed successfully"
    echo "Final memory usage:"
    free -h || echo "free command not available"
    """
}


process methylDackel_extract {
    label 'high_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
    // Add specific error strategy for this process to handle memory issues
    errorStrategy = { task.exitStatus in [143,137,104,134,139] || task.attempt <= 3 ? 'retry' : 'finish' }
    // Increase max retries for this process
    maxRetries = 3
    // Add memory directive to ensure adequate memory allocation
    memory = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def fileSizeGB = md_bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
        // Calculate memory based on file size with a higher multiplier
        def fileBasedMemGB = Math.ceil(fileSizeGB * 2.5).doubleValue() // Increased multiplier from 2.0 to 2.5
        def minMemGB = slurm_profile ? 8.GB : 4.GB // Increased minimum memory

        // Use the maximum of calculated memory or minimum memory, multiplied by attempt number
        def memToUse = Math.max(minMemGB.toGiga(), fileBasedMemGB).GB * task.attempt

        // Cap at max_memory
        check_max(memToUse, 'memory')
    }
    // Add time directive to ensure adequate time allocation
    time = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def default_time = slurm_profile ? 48.h : 12.h // Increased from 24h to 48h for SLURM

        check_max(default_time * task.attempt, 'time')
    }
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"  // Skip procps-ng for macOS
        } else {
            "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8 conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }

    input:
        tuple val(library), path(md_bam), path(md_bai), val(barcodes) 
        path(genome_path)

    output:
        tuple val(library), file('*.methylKit.gz'), emit: extract_output 

    script:
    // Set memory and CPU resources based on BAM file size
    def fileSizeGB = md_bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
    def currentMemoryGB = task.memory.toGiga() // Convert task.memory to GB

    // Ensure all values are explicitly converted to double to avoid type ambiguity
    def maxMemoryGB = params.max_memory.toGiga().doubleValue()
    def currentMemGB = currentMemoryGB.doubleValue()
    def minMemGB = 2.0d
    def fileBasedMemGB = Math.ceil(fileSizeGB * 2.5).doubleValue() // Increased multiplier from 2.0 to 2.5

    // MethylDackel is memory-intensive for large files
    // Scale memory with file size but ensure minimum and respect maximum
    def memoryGB = Math.min(
        maxMemoryGB,
        Math.max(Math.max(currentMemGB, minMemGB), fileBasedMemGB)
    )

    // Adjust CPUs based on available resources and file size
    // Convert to integers for CPU calculations
    def taskCpusInt = task.cpus.intValue()
    def fileSizeInt = Math.ceil(fileSizeGB * 2).intValue()
    def minCpus = 2
    def maxCpus = 8

    def cpusToUse = Math.min(
        taskCpusInt,
        Math.max(minCpus, Math.min(maxCpus, fileSizeInt))
    )

    // Convert cpusToUse to integer for pigz (already an integer from above calculations)
    def cpusToUseInt = cpusToUse

    task.memory = "${memoryGB} GB"

    """
    echo "Input BAM size: ${fileSizeGB} GB"
    echo "Memory allocated for this task: ${task.memory}"
    echo "CPUs allocated: ${cpusToUse}"

    # Monitor memory usage
    echo "Available memory before processing:"
    free -h || echo "free command not available"

    # Create a dedicated temp directory with random name to avoid conflicts
    TEMP_DIR="\${TMPDIR:-${params.tmp_dir}}/methyldackel_\${RANDOM}"
    mkdir -p "\$TEMP_DIR"
    echo "Using temporary directory: \$TEMP_DIR"

    # Check disk space in temp directory
    echo "Disk space in temp directory:"
    df -h "\$TEMP_DIR" || echo "df command not available"

    # Find the genome file
    genome=\$(ls *fa)

    # Step 1: Run MethylDackel extract with memory monitoring
    echo "Step 1: Running MethylDackel extract..."

    # Run MethylDackel with controlled memory usage
    MethylDackel extract --methylKit -q ${params.min_mapq} -@ ${cpusToUseInt} \
        --CHH --CHG -o ${library}.${barcodes} \${genome} ${md_bam}

    # Check exit status of MethylDackel
    methyldackel_exit=\$?
    if [ \$methyldackel_exit -ne 0 ]; then
        echo "MethylDackel extract failed with exit code \$methyldackel_exit"
        echo "This might be due to memory constraints or disk space issues."

        # Check memory usage
        echo "Current memory usage:"
        free -h || echo "free command not available"

        # Check disk space
        echo "Disk space in temp directory:"
        df -h "\$TEMP_DIR" || echo "df command not available"

        # Clean up temp directory
        rm -rf "\$TEMP_DIR"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry with more memory"
            exit 137  # Memory error code that will trigger retry
        fi

        exit \$methyldackel_exit
    fi

    # Verify the output files exist
    if [ ! -f "${library}.${barcodes}_CpG.methylKit" ]; then
        echo "Error: CpG methylKit file is missing"
        exit 1
    fi

    # Step 2: Compress the methylKit files
    echo "Step 2: Compressing methylKit files..."

    # Use pigz for parallel compression
    pigz -p ${cpusToUseInt} *.methylKit

    # Check exit status of pigz
    pigz_exit=\$?
    if [ \$pigz_exit -ne 0 ]; then
        echo "Pigz compression failed with exit code \$pigz_exit"

        # Check disk space
        echo "Disk space in working directory:"
        df -h . || echo "df command not available"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry"
            exit 1
        fi

        exit \$pigz_exit
    fi

    # Verify the compressed files exist
    if [ ! -f "${library}.${barcodes}_CpG.methylKit.gz" ]; then
        echo "Error: Compressed CpG methylKit file is missing"
        exit 1
    fi

    # Clean up temp directory
    rm -rf "\$TEMP_DIR"

    echo "MethylDackel extract process completed successfully"
    echo "Final memory usage:"
    free -h || echo "free command not available"
    """
}
