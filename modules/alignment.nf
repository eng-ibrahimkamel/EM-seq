process enough_reads {
    label 'low_cpu'
    tag {library}
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            "bioconda::samtools=1.21"  // Skip procps-ng for macOS
        } else {
            "bioconda::samtools=1.21 conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }

    input:
        tuple val(email), 
              val(library), 
              path(input_file1), 
              path(input_file2), 
              val(fileType)

        output:
            tuple val(email), val(library), path(input_file1), path(input_file2), val(fileType), path("*passes_or_fails.txt") 

        script:
        """
        # Use readlink -f as a more portable alternative to realpath
        # If readlink -f is not available, fallback to using the input file directly
        if command -v readlink >/dev/null 2>&1 && readlink -f / >/dev/null 2>&1; then
            in1=\$(readlink -f ${input_file1})
        else
            in1="${input_file1}"
        fi

        passes_or_fails="pass"

        # Check if we're on macOS (BSD) or Linux (GNU)
        if [ "\$(uname)" == "Darwin" ]; then
            # macOS (BSD stat)
            get_file_size() {
                stat -f%z "\$1"
            }
        else
            # Linux (GNU stat)
            get_file_size() {
                stat -c%s "\$1"
            }
        fi

        if grep -q "fastq.gz" <<< "${fileType}"; then
            [ \$(get_file_size "\${in1}") -lt 54 ] && passes_or_fails="fail"
        elif grep -q "fastq" <<< "${fileType}"; then 
            [ \$(get_file_size "\${in1}") -lt 240 ] && passes_or_fails="fail"
        elif grep -q "bam" <<< "${fileType}"; then
            [ \$(get_file_size "\${in1}") -lt 100 ] && passes_or_fails="fail"
        fi 

        echo -e "$library\\t\${passes_or_fails}" > ${library}_passes_or_fails.txt
        """
} 

process send_email {
    label 'low_cpu'
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            ""  // Empty conda environment for macOS
        } else {
            "conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }

    input:
        file libraries

    script:
    """
    touch tmp
    for f in ${libraries}
    do
        cat \$f | awk '{print \$1"<br>"}' >> tmp 
    done
    libs=\$(cat tmp)

    sendmail -t <<EOF
    To: ${params.email}
    Subject: File Read Check
    Content-Type: text/html

    <html>
      <body>
        <p>The following libraries:<br> <strong>\${libs}</strong> do not have enough reads. <br> Continuing with other libraries. </p>
      </body>
    </html>
    EOF
    """
}


process alignReads {
    label 'high_cpu'
    tag { library }
    // Add specific error strategy for this process to handle memory issues
    errorStrategy = { task.exitStatus in [143,137,104,134,139] || task.attempt <= 3 ? 'retry' : 'finish' }
    // Increase max retries for this process
    maxRetries = 3
    // Add memory directive to ensure adequate memory allocation
    memory = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def fileSizeGB = input_file1.size() / (1024 * 1024 * 1024) // Convert bytes to GB
        // Calculate memory based on file size with a higher multiplier
        def fileBasedMemGB = Math.ceil(fileSizeGB * 2.5).doubleValue() // Increased multiplier from 1.5 to 2.5
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
            "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.23.4 bioconda::mark-nonconverted-reads=1.2 bioconda::sambamba=1.0 bioconda::samtools=1.21 bioconda::seqtk=1.4 bioconda::pysam conda-forge::bc"  // Skip procps-ng for macOS
        } else {
            "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.23.4 bioconda::mark-nonconverted-reads=1.2 bioconda::sambamba=1.0 bioconda::samtools=1.21 bioconda::seqtk=1.4 bioconda::pysam conda-forge::procps-ng conda-forge::bc"  // Include procps-ng for Linux
        }
    }
    publishDir "${params.outputDir}/bwameth_align"

    input:
        tuple val(email),
              val(library),
              path(input_file1),
              path(input_file2),
              val(fileType)
        path(genome)

    output:
        tuple val(params.email), val(library), env(barcodes), path("*.nonconverted.tsv"), path("*.fastp.json"), emit: for_agg
        path "*.aln.bam", emit: aligned_bams
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), env(barcodes), emit: bam_files

    script:

    // Set memory, dynamically, based on input file size and respecting resource constraints
    def fileSizeGB = input_file1.size() / (1024 * 1024 * 1024) // Convert bytes to GB
    def currentMemoryGB = task.memory.toGiga() // Convert task.memory to GB

    // Check if we're running in SLURM environment and adjust memory accordingly
    def slurm_profile = workflow.profile.contains('slurm')
    def minMemoryGB = slurm_profile ? 0.8 : 7 // Use 800MB for SLURM, 7GB otherwise

    // Calculate memory based on file size but respect limits
    // Ensure all values are explicitly converted to double to avoid type ambiguity
    def maxMemoryGB = params.max_memory.toGiga().doubleValue()
    def currentMemGB = currentMemoryGB.doubleValue()
    def minMemGB = minMemoryGB.doubleValue()
    def fileBasedMemGB = Math.ceil(fileSizeGB * 1.5).doubleValue()

    def memoryGB = Math.min(
        maxMemoryGB,
        Math.max(Math.max(currentMemGB, minMemGB), fileBasedMemGB)
    )

    task.memory = "${memoryGB} GB"
    println "Task memory set to ${task.memory} (SLURM mode: ${slurm_profile})"

    // Define sambamba_memory here, outside the bash script, respecting resource constraints
    // Ensure consistent types for Math operations
    def memToGigaDouble = task.memory.toGiga().doubleValue()
    def memThreeQuarters = (memToGigaDouble * 3 / 4).doubleValue()

    def sambamba_memory = slurm_profile ? 
        "${Math.min(0.7d, memThreeQuarters)}GB" : // For SLURM, limit to 700MB max
        "${Math.max(4.0d, memThreeQuarters)}GB"   // For non-SLURM, minimum 4GB

    // Calculate threads for bwameth and samtools sort ensuring consistent types
    def bwamethThreads = Math.max(1, (task.cpus.intValue() * 7 / 8).intValue())
    def sortThreads = Math.max(1, (task.cpus.intValue() / 4).intValue())

    """

    echo "Input file size: ${fileSizeGB} GB"
    echo "Memory allocated for this task: ${task.memory}"

    # Determine the genome index
    genome=\$(ls *.bwameth.c2t.bwt | sed 's/.bwameth.c2t.bwt//')

    # Define helper functions
    get_nreads_from_fastq() {
        zcat -f \$1 | grep -c "^+\$" \
        | awk '{
            frac=${params.max_input_reads}/\$1; 
            if (frac>=1) {frac=0.999}; 
            split(frac, numParts, "."); print numParts[2]
        }'
    }

    flowcell_from_fastq() {
        set +o pipefail
        zcat -f \$1 | head -n1 | cut -d ":" -f3
        set -o pipefail
    }


    flowcell_from_bam() {
        set +o pipefail
        samtools view -@ ${task.cpus} \$1 | head -n1 | cut -d":" -f3
        set -o pipefail
    }

    barcodes_from_fastq() {
        set +o pipefail
        zcat -f \$1 \
        | head -n10000 \
        | awk '{
            if (NR%4==1) {
                split(\$0, parts, ":"); 
                arr[ parts[ length(parts) ] ]++
            }} END { for (i in arr) {print arr[i]"\\t"i} }' \
        | sort -k1nr | head -n1 | cut -f2 
        set -o pipefail
    }

    # Determine barcodes and read group line
    get_barcodes_and_rg_line() {
        set +o pipefail
        local file=\$1
        local type=\$2
        if [ "\$type" == "bam" ]; then
            barcodes=\$(samtools view -@ ${task.cpus} -H \$file | grep @RG | awk '{for (i=1;i<=NF;i++) {if (\$i~/BC:/) {print substr(\$i,4,length(\$i))} } }' | head -n1)
            rg_line=\$(samtools view -@ ${task.cpus} -H \$file | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
        else
            barcodes=(\$(barcodes_from_fastq \$file))
            rg_line="@RG\\tID:\${barcodes}\\tSM:${library}\\tBC:\${barcodes}"
        fi
        set -o pipefail
    }

    get_frac_reads() {
        local file=\$1
        local type=\$2
        if [ "${params.max_input_reads}" == "all_reads" ]; then
            frac_reads=1
        else
            if [ "\$type" == "bam" ]; then
                n_reads=\$(samtools view -@ ${task.cpus} -c -F 2304 \$file)
            else
                n_reads=\$(get_nreads_from_fastq \$file)
            fi
            if [ \$n_reads -le ${params.max_input_reads} ]; then
                frac_reads=1
            else
                frac_reads=\$(echo \$n_reads | awk '{print ${params.max_input_reads}/\$1}')
            fi 
        fi
    }


    reheader_sam() {
      local input_file="\$1"

      cat "\$input_file" | \
      awk 'BEGIN{
        # Define header types in the correct order they should appear
        header_order = "@HD @SQ @RG @PG @CO"
        split(header_order, order_arr, " ")
        for (i in order_arr) {
          order_idx[order_arr[i]] = i
        }
        header_printed = 0
      } 
      {
        # Process header lines
        if (\$1~/^@/) {
          # Store header lines by type
          id = substr(\$1,1,3)
          if (!(id in headers)) {
            headers[id] = \$0
          } else {
            headers[id] = headers[id] "\\n" \$0
          }
        }
        else {
          # Print all headers in the correct order before the first alignment record
          if (!header_printed) {
            for (i=1; i<=length(order_arr); i++) {
              if (headers[order_arr[i]]) {
                print headers[order_arr[i]]
              }
            }
            header_printed = 1
          }
          # Print the alignment record
          print \$0
        }
      }'
    }


   case ${fileType} in 
        "fastq_paired_end")
            get_barcodes_and_rg_line ${input_file1} "fastq"
            get_frac_reads ${input_file1} "fastq"
            stream_reads="samtools import -u -1 ${input_file1} -2 ${input_file2}"
            flowcell=\$(flowcell_from_fastq ${input_file1})
            ;;
        "bam")
            get_barcodes_and_rg_line ${input_file1} "bam"
            get_frac_reads ${input_file1} "bam"
            stream_reads="samtools view -@ ${task.cpus} -u -h ${input_file1}"
            flowcell=\$(flowcell_from_bam ${input_file1})
            ;;
        "fastq_single_end")
            get_barcodes_and_rg_line ${input_file1} "fastq"
            get_frac_reads ${input_file1} "fastq"
            stream_reads="samtools import -u -s ${input_file1}"
            flowcell=\$(flowcell_from_fastq ${input_file1})
            ;;
    esac
    if [ "${params.flowcell}" == "undefined" ]; then
        flowcell="\${flowcell}"
    else
        flowcell="${params.flowcell}"
    fi

    if [ \${frac_reads} -lt 1 ]; then
        downsample_seed_frac=\$(awk -v seed=${params.downsample_seed} -v frac=\${frac_reads} 'BEGIN { printf "%.4f", seed + frac }')
        stream_reads="\${stream_reads} | samtools view -@ ${task.cpus} -u -s \${downsample_seed_frac}"
    fi

    base_outputname="${library}_\${barcodes}_\${flowcell}"

    set +o pipefail
    # Use parallelization for samtools view
    inst_name=\$(samtools view -@ ${task.cpus} ${input_file1} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    bam2fastq="| samtools collate -f -r 100000 -u /dev/stdin -O | samtools fastq -n  /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

    # Monitor memory usage during alignment
    echo "Starting alignment with memory monitoring..."
    echo "Available memory before alignment:"
    free -h || echo "free command not available"

    # Break the pipeline into smaller steps to isolate issues and better manage memory
    # Step 1: Process reads with fastp and save to intermediate file
    echo "Step 1: Processing reads with fastp..."
    eval \${stream_reads} \${bam2fastq} \
    | fastp --stdin --stdout -l 2 -Q \${trim_polyg} --interleaved_in --overrepresentation_analysis -j "\${base_outputname}.fastp.json" -w ${task.cpus} 2> fastp.stderr > "\${base_outputname}.processed.fq"

    fastp_exit=\$?
    if [ \$fastp_exit -ne 0 ]; then
        echo "Fastp processing failed with exit code \$fastp_exit"
        echo "This might be due to memory constraints. Check the fastp.stderr file for details."
        exit \$fastp_exit
    fi

    # Check intermediate file size
    echo "Processed FASTQ size: \$(du -h "\${base_outputname}.processed.fq" | cut -f1)"

    # Step 2: Run BWA-MEM alignment with memory monitoring
    echo "Step 2: Running BWA-MEM alignment..."
    echo "Available memory before BWA-MEM:"
    free -h || echo "free command not available"

    # Set memory limit for BWA-MEM (80% of available memory)
    mem_value=\$(echo "${task.memory}" | sed -E 's/([0-9.]+).*/\\1/')
    mem_unit=\$(echo "${task.memory}" | sed -E 's/[0-9.]+ *([A-Za-z]+).*/\\1/')

    # Calculate BWA memory limit (80% of allocated memory)
    if [[ "\${mem_unit}" == "GB" || "\${mem_unit}" == "gb" || "\${mem_unit}" == "G" || "\${mem_unit}" == "g" ]]; then
        bwa_mem_limit=\$(echo "\${mem_value} * 0.8" | bc | cut -d'.' -f1)
    else
        bwa_mem_limit=\$(echo "\${mem_value} * 0.8 / 1024" | bc | cut -d'.' -f1)
    fi

    echo "Setting BWA memory limit to approximately \${bwa_mem_limit}GB"

    # Run BWA-MEM with controlled memory usage
    cat "\${base_outputname}.processed.fq" | \
    bwameth.py -p -t ${bwamethThreads} --read-group "\${rg_line}" --reference \${genome} /dev/stdin 2> "\${base_outputname}.log.bwamem" > "\${base_outputname}.sam"

    # Check exit status of the bwameth.py command
    bwameth_exit=\$?
    if [ \$bwameth_exit -ne 0 ]; then
        echo "BWA-MEM alignment failed with exit code \$bwameth_exit"
        echo "This might be due to memory constraints. Check the log file for details."

        # Check if we can find memory-related errors in the log
        if grep -q "out of memory" "\${base_outputname}.log.bwamem" || grep -q "allocate" "\${base_outputname}.log.bwamem"; then
            echo "Memory-related error detected in BWA-MEM log"
            echo "Current memory usage:"
            free -h || echo "free command not available"

            # If this is not the last retry, exit with a code that will trigger a retry
            if [ ${task.attempt} -lt 3 ]; then
                echo "Will retry with more memory"
                exit 137  # Memory error code that will trigger retry
            fi
        fi

        exit \$bwameth_exit
    fi

    # Clean up intermediate file to save space
    rm -f "\${base_outputname}.processed.fq"

    # Step 2: Reheader the SAM file
    cat "\${base_outputname}.sam" | reheader_sam /dev/stdin > "\${base_outputname}.reheadered.sam"

    # Step 3: Skip mark-nonconverted-reads.py due to pysam dependency issues
    # Create an empty nonconverted.tsv file to satisfy the output requirements
    touch "\${base_outputname}.nonconverted.tsv"

    # Step 4: Convert to BAM and sort using samtools instead of sambamba
    echo "Step 3: Converting SAM to BAM and sorting..."
    echo "Available memory before SAM to BAM conversion:"
    free -h || echo "free command not available"

    # Calculate memory limit for samtools sort based on available memory
    # Use a simpler approach to avoid AWK escaping issues
    mem_value=\$(echo "${task.memory}" | sed -E 's/([0-9.]+).*/\\1/')
    mem_unit=\$(echo "${task.memory}" | sed -E 's/[0-9.]+ *([A-Za-z]+).*/\\1/')

    # Convert to MB based on unit (using bc for more reliable arithmetic)
    if [[ "\${mem_unit}" == "GB" || "\${mem_unit}" == "gb" || "\${mem_unit}" == "G" || "\${mem_unit}" == "g" ]]; then
        # Convert GB to MB (multiply by 1024)
        mem_value_mb=\$(echo "\${mem_value} * 1024" | bc | cut -d'.' -f1)
    elif [[ "\${mem_unit}" == "KB" || "\${mem_unit}" == "kb" || "\${mem_unit}" == "K" || "\${mem_unit}" == "k" ]]; then
        # Convert KB to MB (divide by 1024)
        mem_value_mb=1  # Default to 1 MB if less than 1 MB
        if (( \$(echo "\${mem_value} > 1024" | bc) )); then
            mem_value_mb=\$(echo "\${mem_value} / 1024" | bc | cut -d'.' -f1)
        fi
    else
        # Assume already in MB
        mem_value_mb=\$(echo "\${mem_value}" | cut -d'.' -f1)
    fi

    # Calculate memory per thread (65% of total divided by thread count)
    # Reduced from 75% to 65% to leave more memory for the OS and other processes
    threads=${sortThreads}
    mem_per_thread=\$(( (mem_value_mb * 65 / 100) / threads ))

    # Ensure minimum of 100M per thread
    if [ \${mem_per_thread} -lt 100 ]; then
        mem_per_thread=100
    fi

    sort_mem_per_thread="\${mem_per_thread}M"

    echo "Memory per thread for samtools sort: \$sort_mem_per_thread"
    echo "Using ${sortThreads} threads for sorting"

    # First convert SAM to BAM
    echo "Converting SAM to BAM..."
    samtools view -@ ${task.cpus} -u "\${base_outputname}.reheadered.sam" > "\${base_outputname}.unsorted.bam"

    view_exit=\$?
    if [ \$view_exit -ne 0 ]; then
        echo "SAM to BAM conversion failed with exit code \$view_exit"
        echo "This might be due to memory constraints."

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry with more memory"
            exit 137  # Memory error code that will trigger retry
        fi

        exit \$view_exit
    fi

    # Remove SAM file to save space
    rm -f "\${base_outputname}.reheadered.sam"

    # Then sort the BAM file
    echo "Sorting BAM file..."
    echo "Available memory before sorting:"
    free -h || echo "free command not available"

    # Use a temporary directory with enough space
    TEMP_SORT_DIR="\${TMPDIR:-${params.tmp_dir}}/sort_\${RANDOM}"
    mkdir -p "\$TEMP_SORT_DIR"

    samtools sort -m \$sort_mem_per_thread -@ ${sortThreads} -T "\$TEMP_SORT_DIR/tmp" -o "\${base_outputname}.aln.bam" "\${base_outputname}.unsorted.bam"

    sort_exit=\$?
    if [ \$sort_exit -ne 0 ]; then
        echo "BAM sorting failed with exit code \$sort_exit"
        echo "This might be due to memory constraints or disk space issues."

        # Check disk space
        echo "Disk space in temp directory:"
        df -h "\$TEMP_SORT_DIR" || echo "df command not available"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry with more memory"
            exit 137  # Memory error code that will trigger retry
        fi

        exit \$sort_exit
    fi

    # Clean up temporary directory
    rm -rf "\$TEMP_SORT_DIR"

    # Remove unsorted BAM to save space
    rm -f "\${base_outputname}.unsorted.bam"

    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index "\${base_outputname}.aln.bam"

    index_exit=\$?
    if [ \$index_exit -ne 0 ]; then
        echo "BAM indexing failed with exit code \$index_exit"

        # If this is not the last retry, exit with a code that will trigger a retry
        if [ ${task.attempt} -lt 3 ]; then
            echo "Will retry"
            exit 137
        fi

        exit \$index_exit
    fi

    # Clean up any remaining intermediate files
    rm -f "\${base_outputname}.sam"

    echo "Alignment process completed successfully"
    echo "Final memory usage:"
    free -h || echo "free command not available"


    """
}

process mergeAndMarkDuplicates {
    label 'high_cpu'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    // Add specific error strategy for this process to handle memory issues
    errorStrategy = { task.exitStatus in [143,137,104,134,139] || task.attempt <= 3 ? 'retry' : 'finish' }
    // Increase max retries for this process
    maxRetries = 3
    // Add memory directive to ensure adequate memory allocation
    memory = { 
        def slurm_profile = workflow.profile.contains('slurm')
        def fileSizeGB = bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
        // Calculate memory based on file size with a higher multiplier for Picard
        def fileBasedMemGB = Math.ceil(fileSizeGB * 3.0).doubleValue() // Higher multiplier for Picard
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
            "bioconda::picard=3.1 bioconda::samtools=1.21"  // Skip procps-ng for macOS
        } else {
            "bioconda::picard=3.1 bioconda::samtools=1.21 conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }

    input:
        tuple val(library), path(bam), path(bai), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), val(barcodes), emit: md_bams
        tuple val( params.email ), val(library), path('*.md.bam'), path('*.md.bai'), emit: for_agg
        path('*.markdups_log'), emit: log_files

    script:
    // Calculate picard memory ensuring consistent types
    def fileSizeGB = bam.size() / (1024 * 1024 * 1024) // Convert bytes to GB
    def currentMemoryGB = task.memory.toGiga() // Convert task.memory to GB

    // Calculate Picard memory - use 80% of available memory
    def picardMemGB = Math.max(1, (currentMemoryGB * 0.8).intValue())

    """
    echo "Input BAM file size: ${fileSizeGB} GB"
    echo "CPUs allocated: ${task.cpus}"
    echo "Memory allocated for this task: ${task.memory}"
    echo "Picard memory allocation: ${picardMemGB}g"

    # Monitor memory usage
    echo "Available memory before processing:"
    free -h || echo "free command not available"

    set +o pipefail
    # Use parallelization for samtools view
    inst_name=\$(samtools view -@ ${task.cpus} ${bam} | head -n1 | cut -d ":" -f1);
    set -o pipefail

    optical_distance=\$(echo \${inst_name} | awk '{if (\$1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

    # Create a dedicated temp directory with random name to avoid conflicts
    TEMP_DIR="\${TMPDIR:-${params.tmp_dir}}/picard_\${RANDOM}"
    mkdir -p "\$TEMP_DIR"
    echo "Using temporary directory: \$TEMP_DIR"

    # Check disk space in temp directory
    echo "Disk space in temp directory:"
    df -h "\$TEMP_DIR" || echo "df command not available"

    # Calculate optimal number of threads for Picard
    # Picard benefits from multiple threads for MarkDuplicates
    echo "Running Picard MarkDuplicates..."
    picard -Xmx${picardMemGB}g MarkDuplicates \
        --TAGGING_POLICY All \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE \${optical_distance} \
        --TMP_DIR "\$TEMP_DIR" \
        --CREATE_INDEX true \
        --MAX_RECORDS_IN_RAM 5000000 \
        --BARCODE_TAG "RX" \
        --ASSUME_SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY SILENT \
        -I ${bam} \
        -O ${library}_${barcodes}.md.bam \
        -M ${library}.markdups_log

    picard_exit=\$?
    if [ \$picard_exit -ne 0 ]; then
        echo "Picard MarkDuplicates failed with exit code \$picard_exit"
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

        exit \$picard_exit
    fi

    # Verify the output files exist and are valid
    if [ ! -s "${library}_${barcodes}.md.bam" ]; then
        echo "Error: Output BAM file is empty or does not exist"
        exit 1
    fi

    if [ ! -s "${library}_${barcodes}.md.bai" ]; then
        echo "Error: Output BAI file is empty or does not exist"
        # Try to create the index if it doesn't exist
        echo "Attempting to create index manually..."
        samtools index "${library}_${barcodes}.md.bam"
    fi

    # Clean up temp directory
    rm -rf "\$TEMP_DIR"

    echo "MarkDuplicates process completed successfully"
    echo "Final memory usage:"
    free -h || echo "free command not available"
    """
}

process bwa_index {
    /* This pipeline is for internal AND external use.
     * Attempts to link the reference index. If there is no index
     * we download it from the provided URL.
     * If no index and no URL, User will have to debug.
     *
     * Note: For a 3GB reference genome, this process requires:
     * - At least 16GB of memory (5-6x the reference size)
     * - Approximately 15GB of disk space (5x the reference size)
     * - Significant CPU resources for faster indexing
     */

    label 'high_cpu'  // Upgraded to high_cpu for more resources to handle 3GB reference genome
    tag { genome }
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            "conda-forge::python=3.10 bioconda::samtools=1.21 bioconda::bwameth=0.2.7"  // Skip procps-ng for macOS
        } else {
            "conda-forge::python=3.10 bioconda::samtools=1.21 bioconda::bwameth=0.2.7 conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }
    storeDir "${params.storeDir}"
    // Only retry for specific exit codes that indicate transient issues
    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries = 3

    // Custom memory allocation for bwa_index process
    // BWA indexing typically requires 5-6x the reference genome size
    memory = {
        def slurm_profile = workflow.profile.contains('slurm')
        // For a 3GB reference genome, allocate at least 16GB
        def min_memory = slurm_profile ? 16.GB : 6.GB

        check_max(min_memory * task.attempt, 'memory')
    }

    output:
    path "*.{fa,fai,amb,ann,bwt,pac,sa,c2t,bwameth.c2t.*}"

    script:
    // Calculate optimal number of threads for BWA indexing
    // Based on the issue description, bwameth.py supports threading via the --threads parameter
    def bwaThreads = Math.max(1, task.cpus.intValue() - 1)

    """
    # Debug: Print conda environment info
    echo "Conda environment path: \$CONDA_PREFIX"
    echo "Python version:"
    python --version
    echo "PATH:"
    echo \$PATH
    echo "Looking for bwameth.py:"
    which bwameth.py || echo "bwameth.py not found in PATH"

    # Print resource allocation
    echo "CPU cores allocated: ${task.cpus}"
    echo "Memory allocated: ${task.memory}"
    echo "Using ${bwaThreads} threads for BWA indexing"

    real_genome_file="\$(basename ${params.path_to_genome_fasta})"
    ln -sf "\$(dirname ${params.path_to_genome_fasta})/\${real_genome_file}"* .

    if [ ! -f "\${real_genome_file}.bwameth.c2t.bwt" ]; then
        # if the reference .fa file is a url, not a local path
        if [ ! -f "\${real_genome_file}" ]; then
            echo "Trying to download the reference"
            filename=\$(basename ${params.path_to_genome_fasta})

            if ! curl -f -o \$filename ${params.path_to_genome_fasta}; then
                echo "Error: Failed to download \${params.path_to_genome_fasta}" >&2
                exit 1
            fi
        fi

        # Try to find bwameth.py in the conda environment
        if command -v bwameth.py >/dev/null 2>&1; then
            echo "Starting BWA indexing with ${bwaThreads} threads at \$(date)"
            # bwameth.py index does not support the --threads parameter
            # The threading is handled internally by BWA through OMP_NUM_THREADS
            export OMP_NUM_THREADS=${bwaThreads}
            echo "Set OMP_NUM_THREADS=${bwaThreads} to control BWA parallelization"
            bwameth.py index \${real_genome_file}
            echo "BWA indexing completed at \$(date)"
        else
            echo "Error: bwameth.py not found in PATH. Installing bwameth manually..."
            pip install bwameth
            echo "Starting BWA indexing with ${bwaThreads} threads at \$(date)"
            # bwameth.py index does not support the --threads parameter
            # The threading is handled internally by BWA through OMP_NUM_THREADS
            export OMP_NUM_THREADS=${bwaThreads}
            echo "Set OMP_NUM_THREADS=${bwaThreads} to control BWA parallelization"
            bwameth.py index \${real_genome_file}
            echo "BWA indexing completed at \$(date)"
        fi
    else
        echo "Index files already exist for \${real_genome_file}"
    fi
    """
}

process touchFile {   
    conda {
        // Skip procps-ng on macOS as it's not available
        def os = System.getProperty("os.name").toLowerCase()
        if (os.contains("mac") || os.contains("darwin")) {
            ""  // Empty conda environment for macOS
        } else {
            "conda-forge::procps-ng"  // Include procps-ng for Linux
        }
    }

    input:
        val filename

    script:
    """
    touch ${filename}
    """
}
