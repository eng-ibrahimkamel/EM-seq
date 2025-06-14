#!/usr/bin/env bash
# EM-seq Pipeline Test Runner
# This script runs tests for the EM-seq pipeline with different input file types

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
log_info "Starting EM-seq pipeline tests"

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



# Get absolute paths using a more efficient portable_realpath function
get_abs_path() {
    local path="$1"
    # Use built-in functions if available
    if command_exists readlink && readlink -f / >/dev/null 2>&1; then
        readlink -f "$path"
    elif command_exists realpath; then
        realpath "$path"
    else
        # Fallback implementation
        if [[ -d "$path" ]]; then
            (cd "$path" && pwd)
        else
            local dir=$(dirname "$path")
            local base=$(basename "$path")
            echo "$(cd "$dir" && pwd)/$base"
        fi
    fi
}

log_info "Source directory: ${INPUT_DATA_DIR}"
log_info "Destination directory: ${OUTPUT_DIR}"

# Generate test data from minimal set of reads
log_info "Generating test data from minimal set of reads"

mkdir -p "${OUTPUT_DIR}/test"
# Generate fastq files
log_info "Generating uncompressed fastq files"
for read in "R1" "R2"; do
    input_file="${INPUT_DATA_DIR}/emseq-testg_${read}.fastq.gz"
    output_file="${OUTPUT_DIR}/test/emseq-test_${read}.fastq"

    if [ ! -f "${input_file}" ]; then
        log_error "Input file not found: ${input_file}"
        exit 1
    fi

    gunzip -c "${input_file}" > "${output_file}" || {
        log_error "Failed to decompress ${input_file}"
        exit 1
    }
    log_info "Created ${output_file}"
done

# Generate BAM file
log_info "Generating BAM file"
paste -d "\n" <(samtools sort "${OUTPUT_DIR}/test/emseq-test_R1.fastq" | samtools view | awk 'BEGIN{OFS="\t"}{$2=77; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') \
              <(samtools sort "${OUTPUT_DIR}/test/emseq-test_R2.fastq" | samtools view | awk 'BEGIN{OFS="\t"}{$2=141; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') \
| samtools view -u -o "${OUTPUT_DIR}/test/emseq-test.u.bam" || {
    log_error "Failed to generate BAM file"
    exit 1
}
log_info "Created ${OUTPUT_DIR}/test/emseq-test.u.bam"


# List files in the temporary directory
log_info "Files in temporary directory:"
ls -ltr "${OUTPUT_DIR}/test" | while read line; do
    log_info "  ${line}"
done

# Function to run a test with a specific input file pattern
test_pipeline() {
    local file_pattern="$1"
    local test_name="$2"

    log_info "Running test: ${test_name} with pattern: ${file_pattern}"

    # Check if files matching the pattern exist
    local matching_files=$(ls -1 ${file_pattern} 2>/dev/null | wc -l)
    if [ "${matching_files}" -eq 0 ]; then
        log_error "No files found matching pattern: ${file_pattern}"
        return 1
    fi

    log_info "Found ${matching_files} files matching pattern: ${file_pattern}"

    # Determine the correct filename prefix based on test type
    local file_prefix="emseq-testg"
    case "$test_name" in
      bam|fastq)
        file_prefix="emseq-test"
        ;;
    esac

    log_info "Using file prefix: ${file_prefix} for ${test_name} test"

    # Run the Nextflow pipeline with the specified input file
    log_info "Starting Nextflow pipeline for ${test_name}"
    source "${SCRIPT_DIR}/scripts.config"
    pushd ${OUTPUT_DIR}
    nextflow run "${BASE_DIR}/main.nf" \
        --input_glob "${file_pattern}" \
        --path_to_genome_fasta "${GENOME_PATH}" \
        --email "${EMAIL}" \
        --max_input_reads "${MAX_INPUT_READS}" \
        --flowcell "${FLOWCELL}" \
        --storeDir "${STORE_DIR}" \
        --outputDir  "${RESULT_DIR}" \
        -with-report "${METADATA_REPORT_PATH}" \
        -with-timeline "${METADATA_TIMELINE_PATH}" \
        -with-dag "${METADATA_DAG_PATH}" \
        -w "${WORK_DIR}" \
        --read_length ${READ_LENGTH} \
        --enable_neb_agg "${ENABLE_NEB_AGG}" 2>&1 | tee -a "${LOG_FILE}" || {
            log_error "Nextflow pipeline failed for ${test_name}"
            echo "Nextflow pipeline failed for ${test_name}" >> "${LOG_FILE}"
            return 1
        }
    popd

    log_info "Nextflow pipeline completed for ${test_name}"
    echo "Nextflow pipeline succeeded for ${test_name}" >> "${LOG_FILE}"

    # Check results
    log_info "Checking results for ${test_name}"

    # Track check results
    local checks_passed=true

    # Check flagstats
    if cat "${RESULT_DIR}/stats/flagstats/${file_prefix}.flagstat" | grep -q "1972 + 0 properly paired"; then
        log_info "Flagstats check passed for ${test_name}"
        echo "flagstats OK for ${test_name}" >> "${LOG_FILE}"
    else
        log_warning "Flagstats check failed for ${test_name}"
        echo "flagstats not OK for ${test_name}" >> "${LOG_FILE}"
        checks_passed=false
    fi

    # Check alignment metrics
    local alignment_result=$(tail -n2 "${RESULT_DIR}/stats/picard_alignment_metrics/${file_prefix}.alignment_summary_metrics.txt" | \
        awk 'BEGIN{result="alignment metrics not OK"}{if ($1==150 && $3>2200) {result="alignment metrics OK"}}END{print result}')

    if [ "${alignment_result}" = "alignment metrics OK" ]; then
        log_info "Alignment metrics check passed for ${test_name}"
    else
        log_warning "Alignment metrics check failed for ${test_name}"
        checks_passed=false
    fi

    echo "${alignment_result} for ${test_name}" >> "${LOG_FILE}"

    if [ "${checks_passed}" = true ]; then
        return 0
    else
        return 1
    fi
}

# Run tests with different input file types
log_info "Running tests with different input file types"

# Track overall test status
TESTS_PASSED=true

# Test with fastq.gz files
if ! test_pipeline "${INPUT_DATA_DIR}/emseq-test*1.fastq.gz" "fastq_gz"; then
    TESTS_PASSED=false
    log_warning "Test with fastq.gz files failed, but continuing with other tests"
fi

# Test with fastq files
if ! test_pipeline "${INPUT_DATA_DIR}/emseq-test*1.fastq" "fastq"; then
    TESTS_PASSED=false
    log_warning "Test with fastq files failed, but continuing with other tests"
fi

# Test with bam files
if ! test_pipeline "${INPUT_DATA_DIR}/emseq-test*bam" "bam"; then
    TESTS_PASSED=false
    log_warning "Test with bam files failed, but continuing with other tests"
fi


## Display test results
#log_info "Test results:"
#cat "${LOG_FILE}" | while read line; do
#    log_info "  ${line}"
#done

# Clean up the tmp directory (uncomment to enable cleanup)
# log_info "Cleaning up temporary directory"
# if [ -d "${TMP_DIR}" ]; then
#     rm -rf "${TMP_DIR}" || log_warning "Failed to clean up temporary directory"
# fi

# Final status
if [ "${TESTS_PASSED}" = true ]; then
    log_info "All tests completed successfully"
else
    log_warning "Some tests failed, check the log for details"
fi

log_info "Testing completed"
