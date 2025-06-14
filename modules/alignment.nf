process enough_reads {
    label 'low_cpu'
    tag {library}
    conda "bioconda::samtools=1.21"

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
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.23.4 bioconda::mark-nonconverted-reads=1.2 bioconda::sambamba=1.0 bioconda::samtools=1.21 bioconda::seqtk=1.4 bioconda::pysam"
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

    // Set memory, dynamically, based on input file size
    def fileSizeGB = input_file1.size() / (1024 * 1024 * 1024) // Convert bytes to GB
    def currentMemoryGB = task.memory.toGiga() // Convert task.memory to GB
    def memoryGB = Math.max(Math.max(currentMemoryGB, 8), Math.ceil(fileSizeGB * 2)) // Minimum 8GB or 2x file size
    task.memory = "${memoryGB} GB"

    // Define sambamba_memory here, outside the bash script
    def sambamba_memory = "${Math.max(4, (task.memory.toGiga()*3).intdiv(4))}GB"

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
        samtools view \$1 | head -n1 | cut -d":" -f3
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
            barcodes=\$(samtools view -H \$file | grep @RG | awk '{for (i=1;i<=NF;i++) {if (\$i~/BC:/) {print substr(\$i,4,length(\$i))} } }' | head -n1)
            rg_line=\$(samtools view -H \$file | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
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
                n_reads=\$(samtools view -c -F 2304 \$file)
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
            stream_reads="samtools view -u -h ${input_file1}"
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
        stream_reads="\${stream_reads} | samtools view -u -s \${downsample_seed_frac}"
    fi

    base_outputname="${library}_\${barcodes}_\${flowcell}"

    set +o pipefail
    inst_name=\$(samtools view ${input_file1} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    bam2fastq="| samtools collate -f -r 100000 -u /dev/stdin -O | samtools fastq -n  /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

    # Break the pipeline into smaller steps to isolate issues
    # Step 1: Process reads and align
    eval \${stream_reads} \${bam2fastq} \
    | fastp --stdin --stdout -l 2 -Q \${trim_polyg} --interleaved_in --overrepresentation_analysis -j "\${base_outputname}.fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference \${genome} /dev/stdin 2> "\${base_outputname}.log.bwamem" > "\${base_outputname}.sam"

    # Step 2: Reheader the SAM file
    cat "\${base_outputname}.sam" | reheader_sam /dev/stdin > "\${base_outputname}.reheadered.sam"

    # Step 3: Skip mark-nonconverted-reads.py due to pysam dependency issues
    # Create an empty nonconverted.tsv file to satisfy the output requirements
    touch "\${base_outputname}.nonconverted.tsv"

    # Step 4: Convert to BAM and sort using samtools instead of sambamba
    samtools view -u "\${base_outputname}.reheadered.sam" | \
    samtools sort -@ ${Math.max(1,task.cpus.intdiv(4))} -T ${params.tmp_dir}/tmp -o "\${base_outputname}.aln.bam" -

    # Index the BAM file
    samtools index "\${base_outputname}.aln.bam"

    # Clean up intermediate files
    rm -f "\${base_outputname}.sam" "\${base_outputname}.reheadered.sam"


    """
}

process mergeAndMarkDuplicates {
    label 'high_cpu'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    conda "bioconda::picard=3.1 bioconda::samtools=1.21"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), val(barcodes), emit: md_bams
        tuple val( params.email ), val(library), path('*.md.bam'), path('*.md.bai'), emit: for_agg
        path('*.markdups_log'), emit: log_files

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n1 | cut -d ":" -f1);
    set -o pipefail

    optical_distance=\$(echo \${inst_name} | awk '{if (\$1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

    picard -Xmx${task.memory.toGiga()}g MarkDuplicates \
        --TAGGING_POLICY All \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE \${optical_distance} \
        --TMP_DIR ${params.tmp_dir} \
        --CREATE_INDEX true \
        --MAX_RECORDS_IN_RAM 5000000 \
        --BARCODE_TAG "RX" \
        --ASSUME_SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY SILENT \
        -I ${bam} \
        -O ${library}_${barcodes}.md.bam \
        -M ${library}.markdups_log
    """
}

process bwa_index {
    /* This pipeline is for internal AND external use.
     * Attempts to link the reference index. If there is no index
     * we download it from the provided URL.
     * If no index and no URL, User will have to debug.
     */

    label 'low_cpu'
    tag { genome }
    conda "conda-forge::python=3.10 bioconda::samtools=1.21 bioconda::bwameth=0.2.7"
    storeDir "${params.storeDir}"
    errorStrategy = 'retry'
    maxRetries = 3

    output:
    path "*.{fa,fai,amb,ann,bwt,pac,sa,c2t}"

    script:
    """
    # Debug: Print conda environment info
    echo "Conda environment path: \$CONDA_PREFIX"
    echo "Python version:"
    python --version
    echo "PATH:"
    echo \$PATH
    echo "Looking for bwameth.py:"
    which bwameth.py || echo "bwameth.py not found in PATH"

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
            bwameth.py index \${real_genome_file}
        else
            echo "Error: bwameth.py not found in PATH. Installing bwameth manually..."
            pip install bwameth
            bwameth.py index \${real_genome_file}
        fi
    else
        echo "Index files already exist for \${real_genome_file}"
    fi
    """
}

process touchFile {   
    conda "conda-forge::procps-ng"

    input:
        val filename

    script:
    """
    touch ${filename}
    """
}
