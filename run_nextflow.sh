#!/bin/bash
set -e

genome_path="$(pwd)/test_data/reference.fa"

# Use local miniconda installation
export PATH="$(pwd)/bin/miniconda3/bin:$PATH"

# Activate the nextflow environment and run nextflow
source "$(pwd)/bin/miniconda3/bin/activate" nextflow
nextflow main.nf  \
  --input_glob "test_data/*1.fastq.gz" \
  --path_to_genome_fasta "${genome_path}" \
  --email "eng.ibrahimkamel@gmail.com" \
  --max_input_reads 10000 \
  --flowcell "test_pipeline" \
  -with-report  "emseq_metadata_report.html" \
  -with-timeline "emseq_metadata_timeline.html" \
  -with-dag "emseq_metadata_dag.html" \
  -w "$(pwd)/test_data/tmp22" \
  --read_length 151 \
  -resume
