#!/bin/bash

# This script parallelises UMI extraction and trimming across multiple samples
# using a native Bash background worker queue.

# Define directories
INPUT_DIR="d_data"
EXTRACT_DIR="o_outputs/extracted_umis"
OUTPUT_DIR="o_outputs/trimmed_fastqs"

# Create necessary output directories
# mkdir -p "$EXTRACT_DIR"
# mkdir -p "$OUTPUT_DIR"

# Define maximum number of samples to process concurrently
MAX_JOBS=4

# Make an array to hold the sample names for processing
FILE_LIST=(\
          "${INPUT_DIR}/1_beads_1m_lib_S1_R1_001.fastq.gz"\
          # "${INPUT_DIR}/2_beads_2m_lib_S2_R1_001.fastq.gz"\
          # "${INPUT_DIR}/3_beads_3m_lib_S3_R1_001.fastq.gz"\
          # "${INPUT_DIR}/4_beads_1p_lib_S4_R1_001.fastq.gz"\
          # "${INPUT_DIR}/5_beads_2p_lib_S5_R1_001.fastq.gz"\
          # "${INPUT_DIR}/6_beads_3p_lib_S6_R1_001.fastq.gz"\
          # "${INPUT_DIR}/7_input_1m_lib_S7_R1_001.fastq.gz"\
          # "${INPUT_DIR}/8_input_2m_lib_S8_R1_001.fastq.gz"\
          # "${INPUT_DIR}/9_input_3m_lib_S9_R1_001.fastq.gz"\
          # "${INPUT_DIR}/10_input_1p_lib_S10_R1_001.fastq.gz"\
          # "${INPUT_DIR}/11_input_2p_lib_S11_R1_001.fastq.gz"\
          # "${INPUT_DIR}/12_input_3p_lib_S12_R1_001.fastq.gz"\
          # "${INPUT_DIR}/Undetermined_S0_R1_001.fastq.gz"\
          )

# ------------------------------------------------------------------
# Define the core processing pipeline as a background worker function
# ------------------------------------------------------------------
process_sample() {
    local r1_file="$1"
    local r2_file="${r1_file/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    local filename=$(basename "$r1_file")
    local sample_name="${filename/_R1_001.fastq.gz/}"
    
    local umi_r1="${EXTRACT_DIR}/${sample_name}_extracted_R1.fastq.gz"
    local umi_r2="${EXTRACT_DIR}/${sample_name}_extracted_R2.fastq.gz"

    echo "[START] Processing sample: ${sample_name}"

    # Step 1: Extract UMI (8bp) and discard linker (3bp)
    # umi_tools extract \
    #   -I "${r1_file}" \
    #   --read2-in="${r2_file}" \
    #   --bc-pattern=NNNNNNNNXXX \
    #   -S "${umi_r1}" \
    #   --read2-out="${umi_r2}" \
    #   --log="${EXTRACT_DIR}/${sample_name}_umi_extract.log"

    # Step 2: Run fastp on the resulting intermediate files
    # (Using 2 threads per fastp execution as configured)
    fastp \
      --in1 "${umi_r1}" \
      --in2 "${umi_r2}" \
      --out1 "${OUTPUT_DIR}/${sample_name}_trimmed_R1.fastq.gz" \
      --out2 "${OUTPUT_DIR}/${sample_name}_trimmed_R2.fastq.gz" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --length_required 36 \
      --html "${OUTPUT_DIR}/${sample_name}_fastp.html" \
      --json "${OUTPUT_DIR}/${sample_name}_fastp.json" \
      --thread 2 
      # > /dev/null 2>&1 # Silences fastp stdout cluttering the logs

    echo "[DONE] Finished processing sample: ${sample_name}"
}

# ------------------------------------------------------------------
# Main Orchestration Loop
# ------------------------------------------------------------------
echo "Starting parallel preprocessing pipeline (Max concurrent samples: ${MAX_JOBS})..."
echo "--------------------------------------------------"

for r1_file in "${FILE_LIST[@]}"; do
    
    # Invoke the worker function in the background for this sample
    process_sample "$r1_file" &
    
    # Control the queue density: check number of active background tasks
    # 'jobs -p' lists process IDs of background tasks running under this shell
    while [ $(jobs -p | wc -l) -ge "$MAX_JOBS" ]; do
        sleep 10  # Pause for 10 seconds before checking again if a slot has cleared
    done
    
done

# Wait for the very last batch of lingering background tasks to exit cleanly
wait

echo "--------------------------------------------------"
echo "All 13 files parallel-processed successfully!"
echo "Trimmed outputs are in: ${OUTPUT_DIR}"