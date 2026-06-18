#!/bin/bash

# This script aligns your uncompressed, verified FASTQ files 
# directly to the newly generated mouse STAR index.

INDEX_DIR="d_data/refGenome/mouse_star_index"
FASTQ_DIR="o_outputs/trimmed_fastqs"
OUTPUT_DIR="o_outputs/star_alignments"

mkdir -p "$OUTPUT_DIR"

for R1 in "${FASTQ_DIR}"/*_trimmed_R1.fastq
do
    [ -e "$R1" ] || continue

    R1_BASE=$(basename "$R1")
    R2_BASE="${R1_BASE/_trimmed_R1.fastq/_trimmed_R2.fastq}"
    
    R1_PATH="${FASTQ_DIR}/${R1_BASE}"
    R2_PATH="${FASTQ_DIR}/${R2_BASE}"
    
    SAMPLE_PREFIX="${R1_BASE/_trimmed_R1.fastq/}"
    PREFIX="${OUTPUT_DIR}/${SAMPLE_PREFIX}_"
    
    echo "========================================="
    echo "Running direct STAR alignment for: ${SAMPLE_PREFIX}"
    echo "========================================="
    
    # Run STAR with robust header handling flags
    STAR --runThreadN 8 \
         --genomeDir "${INDEX_DIR}" \
         --readFilesType Fastx \
         --readFilesIn "${R1_PATH}" "${R2_PATH}" \
         --outSAMreadID Numbered \
         --outFileNamePrefix "${PREFIX}" \
         --outSAMtype BAM SortedByCoordinate
         
    echo "Finished alignment for ${SAMPLE_PREFIX}."
    echo "-----------------------------------------"
    echo ""
done

echo "All samples successfully aligned!"