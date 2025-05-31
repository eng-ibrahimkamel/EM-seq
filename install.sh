#!/bin/bash
set -e

# Create bin directory if it doesn't exist
mkdir -p bin

# Detect OS type
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    if [[ $(uname -m) == "arm64" ]]; then
        # Apple Silicon (M1/M2)
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    else
        # Intel Mac
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    fi
else
    # Linux and others (default to Linux x86_64)
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
fi

# Download Miniconda3 installer
echo "Downloading Miniconda3 installer..."
curl -O $MINICONDA_URL

# Get the filename from the URL
INSTALLER=$(basename $MINICONDA_URL)

# Make the installer executable
chmod +x $INSTALLER

# Install Miniconda3 to local bin directory (use -b for non-interactive installation)
echo "Installing Miniconda3 to ./bin/miniconda3..."
if [ -d "$(pwd)/bin/miniconda3" ]; then
    echo "Miniconda3 directory already exists. Updating existing installation..."
    ./$INSTALLER -b -u -p $(pwd)/bin/miniconda3
else
    echo "Installing new Miniconda3..."
    ./$INSTALLER -b -p $(pwd)/bin/miniconda3
fi

# Skip conda init to avoid modifying user's shell configuration files
echo "Skipping conda init to avoid modifying shell configuration files..."
# We'll provide instructions for using conda without modifying system files

# Create and activate a new environment for Nextflow
echo "Creating nextflow environment..."
./bin/miniconda3/bin/conda create -n nextflow -y

# We can't directly activate the environment in a script, so we'll use conda run
echo "Installing Nextflow..."
./bin/miniconda3/bin/conda run -n nextflow conda install -c bioconda nextflow -y

# Verify installations
echo "Verifying installations..."
./bin/miniconda3/bin/conda --version
./bin/miniconda3/bin/conda run -n nextflow nextflow -version

# Clean up the installer
echo "Cleaning up..."
rm $INSTALLER

echo "Installation complete!"
echo "To use Nextflow, you can run it directly with:"
echo "  ./bin/miniconda3/bin/conda run -n nextflow nextflow [commands]"
echo ""
echo "Note: We skipped 'conda init' to avoid modifying your shell configuration files."
echo "If you want to use the 'conda activate' command, you would need to either:"
echo "  1. Run './bin/miniconda3/bin/conda init bash' (this will modify your ~/.bash_profile)"
echo "  2. Or use the full path each time: './bin/miniconda3/bin/conda activate nextflow'"
