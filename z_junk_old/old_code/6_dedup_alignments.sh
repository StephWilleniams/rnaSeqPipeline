#!/bin/bash

# This script indexes and deduplicates STAR alignment BAM files using UMI identifiers
# specifically for paired-end RNA-seq datasets (e.g., Takara SMARTer UMI kits).

# Define the directory containing your STAR alignments
ALIGN_DIR="o_outputs/star_alignments"

# Loop directly through all coordinate-sorted BAM files in the directory
for BAM_IN in "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam
do
    # Check if any matching BAM files actually exist to avoid empty glob errors
    [ -e "$BAM_IN" ] || continue

    # Extract the base filename (e.g., sample_1_Aligned.sortedByCoord.out.bam)
    BAM_BASE=$(basename "$BAM_IN")
    
    # Extract the clean sample prefix (e.g., sample_1)
    SAMPLE_PREFIX="${BAM_BASE/_Aligned.sortedByCoord.out.bam/}"
    
    # Define the explicit paths for the output files
    BAM_OUT="${ALIGN_DIR}/${SAMPLE_PREFIX}_deduped.bam"
    STATS_PREFIX="${ALIGN_DIR}/${SAMPLE_PREFIX}_umi_stats"
    
    echo "========================================="
    echo "Processing UMI Deduplication for: ${SAMPLE_PREFIX}"
    echo "Input BAM:  ${BAM_BASE}"
    echo "Output BAM: ${SAMPLE_PREFIX}_deduped.bam"
    echo "========================================="
    
    # 1. Index the incoming STAR BAM file (required by umi_tools)
    echo "Indexing input BAM file..."
    samtools index "$BAM_IN"
    
    # 2. Deduplicate based on UMI + mapping position for paired-end data
    echo "Running umi_tools dedup..."
    umi_tools dedup \
              -I "$BAM_IN" \
              --paired \
              -S "$BAM_OUT" \
              --output-stats="$STATS_PREFIX"
              
    # 3. Index the newly created deduplicated BAM file 
    # (This saves you a step later before running featureCounts)
    echo "Indexing deduplicated BAM file..."
    samtools index "$BAM_OUT"
    
    echo "Finished processing ${SAMPLE_PREFIX}."
    echo ""
done

echo "All BAM files have been successfully indexed and deduplicated!"