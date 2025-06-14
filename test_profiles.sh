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

echo "Testing Nextflow configuration profiles..."

# Test standard profile (local execution)
echo "Testing standard profile (local execution)..."
nextflow config -profile standard | grep "executor = 'local'" && echo "✓ Standard profile correctly sets executor to 'local'" || echo "✗ Standard profile test failed"

# Test SLURM profile
echo "Testing SLURM profile..."
nextflow config -profile slurm | grep "executor = 'slurm'" && echo "✓ SLURM profile correctly sets executor to 'slurm'" || echo "✗ SLURM profile test failed"
nextflow config -profile slurm | grep "queue = 'normal'" && echo "✓ SLURM profile correctly sets queue" || echo "✗ SLURM profile queue setting test failed"
nextflow config -profile slurm | grep "clusterOptions = '--account=your_account'" && echo "✓ SLURM profile correctly sets clusterOptions" || echo "✗ SLURM profile clusterOptions test failed"

echo "Configuration profile tests completed."
