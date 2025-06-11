#!/usr/bin/env bash
set -euo pipefail
set -x

# user migth need custom config file
if [ ! -f nextflow.config ]; then
    echo "Copying example nextflow.config to nextflow.config"
    cp nextflow.config.example nextflow.config
fi

# Run install.sh if Miniconda installation doesn't exist
if [ ! -d "$(pwd)/bin/miniconda3" ]; then
    echo "Running install.sh to set up Miniconda and Nextflow..."
    bash install.sh
fi

# Install samtools in the nextflow environment if needed
if ! ./bin/miniconda3/bin/conda list -n nextflow | grep -q samtools; then
    echo "Installing samtools in the nextflow environment..."
    ./bin/miniconda3/bin/conda install -n nextflow -c conda-forge -c bioconda samtools=1.21 -y
fi

# Use local miniconda installation
export PATH="$(pwd)/bin/miniconda3/bin:$PATH"

# Activate the nextflow environment
source "$(pwd)/bin/miniconda3/bin/activate" nextflow

# Get the current directory
pwd=$(pwd)
tmp="${pwd}/test_data/tmp22"

mkdir -p "${tmp}"

# Set the path to the reference genome
genome_path="${pwd}/test_data/reference.fa"

# Use local miniconda installation
export PATH="${pwd}/bin/miniconda3/bin:$PATH"

# Change to the tmp directory
pushd ${tmp}

# Activate the nextflow environment and run nextflow
source "${pwd}/bin/miniconda3/bin/activate" nextflow

# Run the nextflow workflow
nextflow run ${pwd}/main.nf \
  --input_glob "${pwd}/test_data/*1.fastq.gz" \
  --path_to_genome_fasta "${genome_path}" \
  --email "eng.ibrahimkamel@gmail.com" \
  --max_input_reads 10000 \
  --flowcell "test_pipeline" \
  -with-report "${pwd}/test_data/tmp22/emseq_metadata_report.html" \
  -with-timeline "${pwd}/test_data/tmp22/emseq_metadata_timeline.html" \
  -with-dag "${pwd}/test_data/tmp22/emseq_metadata_dag.html" \
  -w "${pwd}/test_data/tmp22/work" \
  --read_length 151 \
  -resume

# Return to the original directory
popd