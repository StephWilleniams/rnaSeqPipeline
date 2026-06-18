#!/bin/bash
set -e

# --- CONFIGURATION ---
THREADS=8
STAR_INDEX="d_data/refGenome/mouse_star_index/"

# 1. Set your sample number here (e.g., 1, 2, 5, 10)
VAR=1

# 2. Dynamically find the files using the VAR variable and a wildcard (*)
R1_IN=$(ls d_data/${VAR}_*_R1_001.fastq.gz)
R2_IN=$(ls d_data/${VAR}_*_R2_001.fastq.gz)

# 3. Create a dynamic output directory so samples don't overwrite each other!
OUT_DIR="o_outputs/sample_${VAR}"
mkdir -p ${OUT_DIR}

echo "Processing Sample ${VAR}..."
echo "Read 1: ${R1_IN}"
echo "Read 2: ${R2_IN}"
echo "Output going to: ${OUT_DIR}"

# Set the step to skip up to (e.g., SKIP=2 will start running at Step 3)
SKIP=0
# Set the step after which the pipeline should stop (1-7)
STEP=7

echo "========================================================"
echo "Starting Takara SMART-Seq Pipeline (Parallel)"
echo "Configured to skip first ${SKIP} steps and stop after step: ${STEP}"
echo "========================================================"
echo ""

# --- STEP 1: fastp ---
if [[ "$SKIP" -lt 1 ]]; then
    echo "Step 1: Running stupid no fastp umi run..."

    umi_tools extract \
          --extract-method=regex \
          --bc-pattern="(?P<discard_1>.*)(?P<discard_2>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})(?P<discard_3>G{3,5})" \
          -I ${R1_IN} \
          --read2-in=${R2_IN} \
          -S ${OUT_DIR}/extracted_R1.fastq.gz \
          --read2-out=${OUT_DIR}/extracted_R2.fastq.gz \
          --filtered-out=${OUT_DIR}/internal_R1.fastq.gz \
          --filtered-out2=${OUT_DIR}/internal_R2.fastq.gz \
          -L ${OUT_DIR}/umi_tools_extract.log

    if [[ "$STEP" -eq 1 ]]; then echo "Stopping after Step 1."; exit 0; fi
else
    echo "Skipping Step 1..."
fi

# --- STEP 2: Official UMI Extraction ---
if [[ "$SKIP" -lt 2 ]]; then
    echo "Step 2: fastping extracting UMIs..."

    fastp -i ${OUT_DIR}/extracted_R1.fastq.gz -I ${OUT_DIR}/extracted_R2.fastq.gz \
          -o ${OUT_DIR}/extracted-2_R1.fastq.gz -O ${OUT_DIR}/extracted-2_R2.fastq.gz \
          --html ${OUT_DIR}/fastp_report.html --json ${OUT_DIR}/fastp_report.json \
          --thread ${THREADS} 2>${OUT_DIR}/fastp_run.log

    if [[ "$STEP" -eq 2 ]]; then echo "Stopping after Step 2."; exit 0; fi
else
    echo "Skipping Step 2..."
fi

# --- STEP 3: Parallel Alignment (STAR) ---
if [[ "$SKIP" -lt 3 ]]; then
    echo "Step 3a: Aligning 5' UMI reads..."
    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${OUT_DIR}/extracted-2_R1.fastq.gz ${OUT_DIR}/extracted-2_R2.fastq.gz \
         --readFilesCommand "gzip -d -c" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${OUT_DIR}/umi_

    # echo "Step 3b: Aligning Internal reads..."
    # STAR --runThreadN ${THREADS} \
    #      --genomeDir ${STAR_INDEX} \
    #      --readFilesIn ${OUT_DIR}/internal_R1.fastq.gz ${OUT_DIR}/internal_R2.fastq.gz \
    #      --readFilesCommand "gzip -d -c" \
    #      --outSAMtype BAM SortedByCoordinate \
    #      --outFileNamePrefix ${OUT_DIR}/internal_

    if [[ "$STEP" -eq 3 ]]; then echo "Stopping after Step 3."; exit 0; fi
else
    echo "Skipping Step 3..."
fi

# --- STEP 4: Indexing Initial BAMs ---
if [[ "$SKIP" -lt 4 ]]; then
    echo "Step 4: Indexing BAM files..."
    samtools index ${OUT_DIR}/umi_Aligned.sortedByCoord.out.bam
    # samtools index ${OUT_DIR}/internal_Aligned.sortedByCoord.out.bam

    if [[ "$STEP" -eq 4 ]]; then echo "Stopping after Step 4."; exit 0; fi
else
    echo "Skipping Step 4..."
fi

# --- STEP 5: Parallel Deduplication ---
if [[ "$SKIP" -lt 5 ]]; then
    echo "Step 5a: Deduplicating 5' UMI BAM with umi_tools..."
    umi_tools dedup -I ${OUT_DIR}/umi_Aligned.sortedByCoord.out.bam \
                    --output-stats=${OUT_DIR}/umi_dedup_stats \
                    --paired \
                    --chimeric-pairs=discard \
                    --unpaired-reads=discard \
                    -S ${OUT_DIR}/umi_deduplicated.bam \
                    -L ${OUT_DIR}/umi_tools_extract.log

    # echo "Step 5b: Deduplicating Internal BAM with samtools..."
    # samtools collate -O -u -@ ${THREADS} ${OUT_DIR}/internal_Aligned.sortedByCoord.out.bam | \
    # samtools fixmate -m -u - - | \
    # samtools sort -u -@ ${THREADS} | \
    # samtools markdup -r -@ ${THREADS} - ${OUT_DIR}/internal_deduplicated.bam

    if [[ "$STEP" -eq 5 ]]; then echo "Stopping after Step 5."; exit 0; fi
else
    echo "Skipping Step 5..."
fi

# --- STEP 6: Final Production Indexing ---
if [[ "$SKIP" -lt 6 ]]; then
    echo "Step 6: Indexing final deduplicated BAM files..."
    samtools index ${OUT_DIR}/umi_deduplicated.bam
    # samtools index ${OUT_DIR}/internal_deduplicated.bam

    if [[ "$STEP" -eq 6 ]]; then echo "Stopping after Step 6."; exit 0; fi
else
    echo "Skipping Step 6..."
fi

# --- STEP 7: Gene Counting (featureCounts) ---
if [[ "$SKIP" -lt 7 ]]; then
    echo "Step 7: Counting reads per gene with featureCounts..."
    featureCounts -p -T ${THREADS} \
                  -a d_data/refGenome/genomic.gtf \
                  -t exon \
                  -g gene_id \
                  -o ${OUT_DIR}/gene_counts.txt \
                  ${OUT_DIR}/umi_deduplicated.bam
                  #{OUT_DIR}/umi_deduplicated.bam ${OUT_DIR}/internal_deduplicated.bam

    if [[ "$STEP" -eq 7 ]]; then echo "Stopping after Step 7."; exit 0; fi
else
    echo "Skipping Step 7..."
fi

echo "========================================================"
echo "Pipeline complete! Final counts matrix generated."
echo "========================================================"