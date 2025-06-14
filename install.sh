#!/bin/bash
# EM-seq Pipeline Installer
# This script installs Miniconda and Nextflow for the EM-seq pipeline

# Enable error handling
set -e
set -o pipefail

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

# Start installation
log_info "Starting EM-seq pipeline installation"

# Detect OS type and architecture
log_info "Detecting operating system and architecture"
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    if [[ $(uname -m) == "arm64" ]]; then
        # Apple Silicon (M1/M2)
        MINICONDA_URL="${MINICONDA_URL_MACOS_ARM64}"
        log_info "Detected macOS on Apple Silicon (ARM64)"
    else
        # Intel Mac
        MINICONDA_URL="${MINICONDA_URL_MACOS_X86_64}"
        log_info "Detected macOS on Intel (x86_64)"
    fi
else
    # Linux and others
    if [[ $(uname -m) == "aarch64" || $(uname -m) == "arm64" ]]; then
        # ARM64 architecture
        MINICONDA_URL="${MINICONDA_URL_LINUX_ARM64}"
        log_info "Detected Linux on ARM64 architecture"
    else
        # Default to x86_64
        MINICONDA_URL="${MINICONDA_URL_LINUX_X86_64}"
        log_info "Detected Linux on x86_64 architecture"
    fi
fi

# Download Miniconda3 installer with progress bar
log_info "Downloading Miniconda3 installer from ${MINICONDA_URL}"
if command_exists curl; then
    curl -# -L -o ${BIN_SETUP_DIR}/miniconda.sh "${MINICONDA_URL}" || {
        log_error "Failed to download Miniconda installer"
        exit 1
    }
    INSTALLER="${BIN_SETUP_DIR}/miniconda.sh"
else
    log_error "curl command not found. Please install curl and try again."
    exit 1
fi

# Verify the downloaded file
if ! file_exists "${INSTALLER}"; then
    log_error "Installer file not found after download"
    exit 1
fi

# Make the installer executable
log_info "Making installer executable"
chmod +x "${INSTALLER}"

# Install Miniconda3 to local bin directory (use -b for non-interactive installation)
if dir_exists "${MINICONDA_SETUP_DIR}"; then
    log_info "Miniconda3 directory already exists. Updating existing installation..."
    "${INSTALLER}" -b -u -p "${MINICONDA_SETUP_DIR}" || {
        log_error "Failed to update Miniconda"
        exit 1
    }
else
    log_info "Installing new Miniconda3 to ${MINICONDA_SETUP_DIR}..."
    "${INSTALLER}" -b -p "${MINICONDA_SETUP_DIR}" || {
        log_error "Failed to install Miniconda"
        exit 1
    }
fi

# Skip conda init to avoid modifying user's shell configuration files
log_info "Skipping conda init to avoid modifying shell configuration files"

# Create a new environment for Nextflow
log_info "Creating ${CONDA_ENV_NAME} environment"
"${MINICONDA_SETUP_DIR}/bin/conda" create -n "${CONDA_ENV_NAME}" -y || {
    log_error "Failed to create conda environment"
    exit 1
}

# Install Nextflow in the environment
log_info "Installing Nextflow and required dependencies in the ${CONDA_ENV_NAME} environment"
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS - skip procps-ng which is not available for macOS
    "${MINICONDA_SETUP_DIR}/bin/conda" run -n "${CONDA_ENV_NAME}" conda install -c conda-forge -c bioconda nextflow samtools=1.21 -y || {
        log_error "Failed to install Nextflow and dependencies"
        exit 1
    }
else
    # Linux - include procps-ng
    "${MINICONDA_SETUP_DIR}/bin/conda" run -n "${CONDA_ENV_NAME}" conda install -c conda-forge -c bioconda nextflow samtools=1.21 procps-ng -y || {
        log_error "Failed to install Nextflow and dependencies"
        exit 1
    }
}

# Verify installations
log_info "Verifying installations"
CONDA_VERSION=$("${MINICONDA_SETUP_DIR}/bin/conda" --version)
log_info "Conda version: ${CONDA_VERSION}"

NEXTFLOW_VERSION=$("${MINICONDA_SETUP_DIR}/bin/conda" run -n "${CONDA_ENV_NAME}" nextflow -version | head -n 1)
log_info "Nextflow version: ${NEXTFLOW_VERSION}"

# Clean up the installer
log_info "Cleaning up temporary files"
rm "${INSTALLER}"

# Installation complete
log_info "Installation complete!"
log_info "To use Nextflow, you can run it directly with:"
log_info "  ${MINICONDA_SETUP_DIR}/bin/conda run -n ${CONDA_ENV_NAME} nextflow [commands]"
log_info ""
log_info "Note: We skipped 'conda init' to avoid modifying your shell configuration files."
log_info "If you want to use the 'conda activate' command, you would need to either:"
log_info "  1. Run '${MINICONDA_SETUP_DIR}/bin/conda init bash' (this will modify your ~/.bash_profile)"
log_info "  2. Or use the full path each time: '${MINICONDA_SETUP_DIR}/bin/conda activate ${CONDA_ENV_NAME}'"
