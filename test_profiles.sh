#!/bin/bash
set -euo pipefail

echo "Testing Nextflow configuration profiles..."

# Test standard profile (local execution)
echo "Testing standard profile (local execution)..."
nextflow config -profile standard | grep "executor = 'local'" && echo "✓ Standard profile correctly sets executor to 'local'" || echo "✗ Standard profile test failed"

# Test SLURM profile
echo "Testing SLURM profile..."
nextflow config -profile slurm | grep "executor = 'slurm'" && echo "✓ SLURM profile correctly sets executor to 'slurm'" || echo "✗ SLURM profile test failed"
nextflow config -profile slurm | grep "queue = 'general'" && echo "✓ SLURM profile correctly sets queue" || echo "✗ SLURM profile queue setting test failed"
nextflow config -profile slurm | grep "clusterOptions = '--account=your_account'" && echo "✓ SLURM profile correctly sets clusterOptions" || echo "✗ SLURM profile clusterOptions test failed"

echo "Configuration profile tests completed."
