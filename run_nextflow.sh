#!/usr/bin/env bash
# EM-seq Pipeline Runner
# This script runs the Nextflow pipeline for EM-seq data analysis

# Enable error handling
set -euo pipefail

# Source the configuration file
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts.config"

# Function to handle errors
handle_error() {
    local exit_code=$?
    local line_number=$1
    log_error "Error on line ${line_number}: Command exited with status ${exit_code}"
    exit ${exit_code}
}

# Set up error trap
trap 'handle_error $LINENO' ERR

# Start execution
log_info "Starting EM-seq Nextflow pipeline execution"

# Check for nextflow.config file
if [ ! -f "${BASE_DIR}/nextflow.config" ]; then
    log_info "Nextflow config file not found, creating from example"
    if [ -f "${BASE_DIR}/nextflow.config.example" ]; then
        cp "${BASE_DIR}/nextflow.config.example" "${BASE_DIR}/nextflow.config"
        log_info "Created nextflow.config from example"
    else
        log_error "nextflow.config.example not found. Cannot create config file."
        exit 1
    fi
fi

# Check for Miniconda installation
if [ ! -d "${MINICONDA_SETUP_DIR}" ]; then
    log_info "Miniconda not found, running installer"
    bash "${BASE_DIR}/install.sh"
    if [ ! -d "${MINICONDA_SETUP_DIR}" ]; then
        log_error "Installation failed: Miniconda directory not found after running install.sh"
        exit 1
    fi
    log_info "Miniconda installation completed"
fi

# Use local miniconda installation
export PATH="${MINICONDA_SETUP_DIR}/bin:$PATH"

# Activate the nextflow environment
log_info "Activating Nextflow environment"
source "${MINICONDA_SETUP_DIR}/bin/activate" "${CONDA_ENV_NAME}"



# Check if genome file exists
if ! file_exists "${GENOME_PATH}"; then
    log_error "Genome file not found: ${GENOME_PATH}"
    exit 1
fi

# Check if input files exist
INPUT_FILES_COUNT=$(ls -1 ${INPUT_GLOB_PATTERN} 2>/dev/null | wc -l)
if [ "${INPUT_FILES_COUNT}" -eq 0 ]; then
    log_warning "No input files found matching pattern: ${INPUT_GLOB_PATTERN}"
    log_warning "Pipeline may fail if no input files are available"
else
    log_info "Found ${INPUT_FILES_COUNT} input files matching pattern: ${INPUT_GLOB_PATTERN}"
fi

## Change to the tmp directory
#log_info "Changing to temporary directory"
#pushd "${TMP_DIR_NEXTFLOW}" || {
#    log_error "Failed to change to directory: ${TMP_DIR_NEXTFLOW}"
#    exit 1
#}

# Run the nextflow workflow
log_info "Running Nextflow workflow"
log_info "Command: nextflow run ${BASE_DIR}/main.nf with parameters:"
log_info "  --input_glob: ${INPUT_GLOB_PATTERN}"
log_info "  --path_to_genome_fasta: ${GENOME_PATH}"
log_info "  --email: ${EMAIL}"
log_info "  --max_input_reads: ${MAX_INPUT_READS}"
log_info "  --flowcell: ${FLOWCELL}"
log_info "  --read_length: ${READ_LENGTH}"

pushd ${OUTPUT_DIR}
nextflow run "${BASE_DIR}/main.nf" \
  --input_glob "${INPUT_GLOB_PATTERN}" \
  --path_to_genome_fasta "${GENOME_PATH}" \
  --email "${EMAIL}" \
  --max_input_reads ${MAX_INPUT_READS} \
  --flowcell "${FLOWCELL}" \
  --storeDir "${STORE_DIR}" \
  -with-report "${METADATA_REPORT_PATH}" \
  -with-timeline "${METADATA_TIMELINE_PATH}" \
  --outputDir  "${RESULT_DIR}" \
  -with-dag "${METADATA_DAG_PATH}" \
  -w "${WORK_DIR}" \
  --read_length ${READ_LENGTH} \
  -resume || {
    log_error "Nextflow pipeline execution failed"
    popd || true
    exit 1
  }
popd
log_info "Nextflow pipeline execution completed successfully"


log_info "Pipeline execution finished. Results are available in: ${RESULT_DIR}"
log_info "Report: ${METADATA_REPORT_PATH}"
log_info "Timeline: ${METADATA_TIMELINE_PATH}"
log_info "DAG: ${METADATA_DAG_PATH}"
