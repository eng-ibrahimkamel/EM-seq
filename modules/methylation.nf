
process methylDackel_mbias {
    label 'medium_cpu'
    errorStrategy 'retry'
    tag "${library}"
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
    def fileBasedMemGB = Math.ceil(fileSizeGB * 1.2).doubleValue()

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
    genome=\$(ls *fa)
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ${library}_${barcodes}_combined_mbias.tsv
    # Use parallelization for samtools view
    chrs=(`samtools view -@ ${cpusToUse} -H ${md_bam} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

    for chr in \${chrs[*]}; do
        for context in CHH CHG CpG; do
            arg=''
            if [ "\$context" = 'CHH' ]; then
            arg='--CHH --noCpG'
            elif [ "\$context" = 'CHG' ]; then
            arg='--CHG --noCpG'
            fi
            # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
            # not sure why we need both --keepDupes and -F, probably a bug in mbias
            join -t \$'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                MethylDackel mbias --noSVG \$arg -@ ${cpusToUse} -r \$chr \${genome} ${md_bam} | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            <( \
                MethylDackel mbias --noSVG --keepDupes -F 2816 \$arg -@ ${cpusToUse} -r \$chr \${genome} ${md_bam} | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            | sed "s/^/\${chr}\t\${context}\t/" \
            >> ${library}_${barcodes}_combined_mbias.tsv
        done
    done
    # makes the svg files for trimming checks
    MethylDackel mbias -@ ${cpusToUse} --noCpG --CHH --CHG -r \${chrs[0]} \${genome} ${md_bam} ${library}_chn
    # Check OS type and use appropriate sed syntax
    if [[ "\$(uname)" == "Darwin" ]]; then
        # macOS
        for f in *chn*.svg; do sed -i '' "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f; done;
    else
        # Linux and other Unix-like systems
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f; done;
    fi

    MethylDackel mbias -@ ${cpusToUse} -r \${chrs[0]} \${genome} ${md_bam} ${library}_cpg
    # Check OS type and use appropriate sed syntax
    if [[ "\$(uname)" == "Darwin" ]]; then
        # macOS
        for f in *cpg*.svg; do sed -i '' "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f; done;
    else
        # Linux and other Unix-like systems
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f; done;
    fi
    """
}


process methylDackel_extract {
    label 'high_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
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
    def fileBasedMemGB = Math.ceil(fileSizeGB * 2).doubleValue()

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

    genome=\$(ls *fa)
    MethylDackel extract --methylKit -q 20 -@ ${cpusToUseInt} \
        --CHH --CHG -o ${library}.${barcodes} \${genome} ${md_bam} 
    pigz -p ${cpusToUseInt} *.methylKit 
    """
}
