#!/bin/bash

### Title: RNA-Seq Pipeline 
### Authors: Emily Skates, Stephen Williams
## This script implements a full RNA-Seq processing workflow for Takara SMART-Seq data, including trimming, UMI extraction, alignment, deduplication, and gene (feature) counting. 
## The file `bioqc_environment.yml` should be used to set up the conda environment with all necessary tools (umi_tools, fastp, STAR, samtools, featureCounts, etc.).

## Note: sections of this script are commented out, such as the internal read processing, to focus on the primary 5' UMI read workflow. 
## For now, these can be uncommented and adapted as needed for additional processing of internal reads.
## At some point I may rewrite the code to be a bit more flexible and have this as an toggle.

# --- ARGUMENT CHECK ---
# Check if a sample number was provided when the script was run
if [ -z "$1" ]; then
    echo "ERROR: No sample number provided."
    echo "Usage: ./code.sh <sample_number>"
    echo "Example: ./code.sh 5"
    exit 1
fi

# --- CONFIGURATION ---
set -e # Exit immediately if a command exits with a non-zero status
THREADS=8 # Number of threads to use for parallel processing
STAR_INDEX="d_data/refGenome/mouse_star_index/" # Path to STAR genome index directory

# Set VAR to the first command-line argument
VAR=$1 

# --- Directory Setup ---
R1_IN=$(ls d_data/${VAR}_*_R1_001.fastq.gz) # Input FASTQ for Read 1
R2_IN=$(ls d_data/${VAR}_*_R2_001.fastq.gz) # Input FASTQ for Read 2
OUT_DIR="o_outputs/sample_${VAR}" # Output directory for processed files
mkdir -p ${OUT_DIR} # Create output directory if it doesn't exist

# Set the step to skip to (e.g., SKIP=2 starts at Step 3) and the step to stop after (1-7)
SKIP=0
STEP=7

echo ""
echo "========================================================"
echo "Starting Pipeline"
echo "Processing Sample ${VAR}..."
echo "Read 1: ${R1_IN}"
echo "Read 2: ${R2_IN}"
echo "Output going to: ${OUT_DIR}"
echo "========================================================"
echo ""

# --- STEP 1: Extract UMIs ---

# INPUTS:
# --extract-method=regex: Use regex-based UMI extraction
# --bc-pattern: Define the regex pattern for UMI extraction
#   - (?P<discard_1>.*): Discard any leading sequence (non-greedy)
#   - (?P<discard_2>ATTGCGCAATG){s<=2}: Match the adapter sequence with up to 2 mismatches
#   - (?P<umi_1>.{8}): Capture the 8bp UMI sequence
#   - (?P<discard_3>G{5}|G{3}): Discard the trailing poly-G sequence (5 or 3 Gs)
# -I: Input FASTQ for Read 1
# --read2-in: Input FASTQ for Read 2
# -S: Output FASTQ for extracted Read 1
# --read2-out: Output FASTQ for extracted Read 2
# --filtered-out: Output FASTQ for discarded Read 1 (internal reads)
# --filtered-out2: Output FASTQ for discarded Read 2 (internal reads)
# -L: Log file for umi_tools extract

# OUTPUTS:
# - ${OUT_DIR}/extracted_R1.fastq.gz: FASTQ with extracted UMIs for Read 1
# - ${OUT_DIR}/extracted_R2.fastq.gz: FASTQ with extracted UMIs for Read 2
# - ${OUT_DIR}/internal_R1.fastq.gz: FASTQ with discarded internal reads for Read 1
# - ${OUT_DIR}/internal_R2.fastq.gz: FASTQ with discarded internal reads for Read 2
# - ${OUT_DIR}/umi_tools_extract.log: Log file for umi_tools extract

# if [[ "$SKIP" -lt 1 ]]; then
#     echo "Step 1: Running umi extraction..."
#     echo ""
#     umi_tools extract \
#           --extract-method=regex \
#           --bc-pattern="(?P<discard_1>.{0,10})(?P<discard_2>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})(?P<discard_3>G{5}|G{3})" \
#           -I ${R1_IN} \
#           --read2-in=${R2_IN} \
#           -S ${OUT_DIR}/extracted_R1.fastq.gz \
#           --read2-out=${OUT_DIR}/extracted_R2.fastq.gz \
#           --filtered-out=${OUT_DIR}/internal_R1.fastq.gz \
#           --filtered-out2=${OUT_DIR}/internal_R2.fastq.gz \
#           -L ${OUT_DIR}/umi_tools_extract.log
#     if [[ "$STEP" -eq 1 ]]; then echo "Stopping after Step 1."; exit 0; fi
# else
#     echo "Skipping Step 1..."
#     echo ""
# fi

# --- STEP 1: Extract UMIs ---
# Note: this is the new C++ version of the UMI extraction step, which replaces the original umi_tools command. 

# INPUTS:
# -r: Input FASTQ for Read 1
# -R: Input FASTQ for Read 2
# -o: Output FASTQ for extracted Read 1
# -O: Output FASTQ for extracted Read 2
# -L: Log file for umi_tools extract
# -i: Output FASTQ for discarded internal reads for Read 1
# -I: Output FASTQ for discarded internal reads for Read 2

# OUTPUTS:
# - ${OUT_DIR}/extracted_R1.fastq.gz: FASTQ with extracted UMIs for Read 1
# - ${OUT_DIR}/extracted_R2.fastq.gz: FASTQ with extracted UMIs for Read 2
# - ${OUT_DIR}/internal_R1.fastq.gz: FASTQ with discarded internal reads for Read 1
# - ${OUT_DIR}/internal_R2.fastq.gz: FASTQ with discarded internal reads for Read 2
# - ${OUT_DIR}/cumi_progress.log: Log file for umi extraction progress and statistics

if [[ "$SKIP" -lt 1 ]]; then
    echo "Step 1: Running umi extraction..."
    echo ""
    ./c_code/cumi_tools -r ${R1_IN} \
                 -R ${R2_IN} \
                 -o ${OUT_DIR}/extracted_R1.fastq.gz \
                 -O ${OUT_DIR}/extracted_R2.fastq.gz \
                 -L ${OUT_DIR}/cumi_progress.log \
                 -i ${OUT_DIR}/internal_R1.fastq.gz \
                 -I ${OUT_DIR}/internal_R2.fastq.gz
    if [[ "$STEP" -eq 1 ]]; then echo "Stopping after Step 1."; exit 0; fi
else
    echo "Skipping Step 1..."
    echo ""
fi

# --- STEP 2: Fastp trimming ---

# INPUTS:
# -i: Input FASTQ for Read 1 (extracted UMIs)
# -I: Input FASTQ for Read 2 (extracted UMIs)
# -o: Output FASTQ for trimmed Read 1
# -O: Output FASTQ for trimmed Read 2
# --html: Output HTML report for fastp
# --json: Output JSON report for fastp
# --thread: Number of threads to use for fastp

# OUTPUTS:
# - ${OUT_DIR}/extracted-2_R1.fastq.gz: Trimmed FASTQ for Read 1
# - ${OUT_DIR}/extracted-2_R2.fastq.gz: Trimmed FASTQ for Read 2
# - ${OUT_DIR}/fastp_report.html: HTML report generated by fastp
# - ${OUT_DIR}/fastp_report.json: JSON report generated by fastp
# - ${OUT_DIR}/fastp_run.log: Log file for fastp execution

if [[ "$SKIP" -lt 2 ]]; then
    echo "Step 2: fastp trimming extracted sequences..."
    echo ""
    fastp -i ${OUT_DIR}/extracted_R1.fastq.gz -I ${OUT_DIR}/extracted_R2.fastq.gz \
          -o ${OUT_DIR}/extracted-2_R1.fastq.gz -O ${OUT_DIR}/extracted-2_R2.fastq.gz \
          --html ${OUT_DIR}/fastp_report.html --json ${OUT_DIR}/fastp_report.json \
          --thread ${THREADS} 2>${OUT_DIR}/fastp_run.log
    if [[ "$STEP" -eq 2 ]]; then echo "Stopping after Step 2."; exit 0; fi
else
    echo "Skipping Step 2..."
    echo ""
fi

# --- STEP 3: Alignment (STAR) ---

# Note: Star alignment can be problematic, particularly the release version available through bioconda.
# To overcome this issue we use version 2.10.7b, which is not available through normal bioconda.
# To get it run the following command:
# CONDA_SUBDIR=osx-64 mamba install -y bioconda::star=2.7.10b
# This version can then cause environment issues, which prevent the use of other tools. Run:
# conda clean --all -y
# conda remove --force star -y
# Temporarily handle other environment issue and then re-install star if this problem occurs.

# INPUTS:
# --runThreadN: Number of threads to use for STAR alignment
# --genomeDir: Path to STAR genome index directory
# --readFilesIn: Input FASTQ files for Read 1 and Read 2 (trimmed with fastp)
# --readFilesCommand: Command to decompress input FASTQ files (gzip -d -c)
# --outSAMtype: Output format for aligned reads (BAM sorted by coordinate)
# --outFileNamePrefix: Prefix for output files generated by STAR

# OUTPUTS:
# - ${OUT_DIR}/star_Aligned.sortedByCoord.out.bam: BAM file with aligned reads sorted by coordinate
# - ${OUT_DIR}/star_Log.out: Log file with STAR alignment summary and statistics
# - ${OUT_DIR}/star_Log.final.out: Log file with detailed STAR alignment results and metrics

if [[ "$SKIP" -lt 3 ]]; then
    echo "Step 3a: Aligning 5' UMI reads..."
    echo ""
    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${OUT_DIR}/extracted-2_R1.fastq.gz ${OUT_DIR}/extracted-2_R2.fastq.gz \
         --readFilesCommand "gzip -d -c" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${OUT_DIR}/umi_

    # --- INTERNAL READS ---
    # echo "Step 3b: Aligning Internal reads..."
    # STAR --runThreadN ${THREADS} \
    #      --genomeDir ${STAR_INDEX} \
    #      --readFilesIn ${OUT_DIR}/internal_R1.fastq.gz ${OUT_DIR}/internal_R2.fastq.gz \
    #      --readFilesCommand "gzip -d -c" \
    #      --outSAMtype BAM SortedByCoordinate \
    #      --outFileNamePrefix ${OUT_DIR}/internal_
    # -----------------------

    if [[ "$STEP" -eq 3 ]]; then echo "Stopping after Step 3."; exit 0; fi
else
    echo "Skipping Step 3..."
    echo ""
fi

# --- STEP 4: Indexing Initial BAMs ---

# INPUTS:
# -I: Input BAM file with aligned reads sorted by coordinate

# OUTPUTS:
# - ${OUT_DIR}/umi_Aligned.sortedByCoord.out.bam.bai: BAM index file for the aligned reads

if [[ "$SKIP" -lt 4 ]]; then
    echo "Step 4: Indexing BAM files..."
    echo ""
    samtools index ${OUT_DIR}/umi_Aligned.sortedByCoord.out.bam

    # --- INTERNAL READS ---
    # samtools index ${OUT_DIR}/internal_Aligned.sortedByCoord.out.bam
    # ----------------------

    if [[ "$STEP" -eq 4 ]]; then echo "Stopping after Step 4."; exit 0; fi
else
    echo "Skipping Step 4..."
    echo ""
fi

# --- STEP 5: Deduplication ---

# Possible change: https://peerj.com/articles/8275/

# INPUTS:
# -I: Input BAM file with aligned reads sorted by coordinate
# --output-stats: Output file for deduplication statistics
# --paired: Indicate that the input BAM contains paired-end reads
# --chimeric-pairs=discard: Discard chimeric read pairs during deduplication
# --unpaired-reads=discard: Discard unpaired reads during deduplication
# -S: Output BAM file with deduplicated reads
# -L: Log file for umi_tools deduplication

# OUTPUTS:
# - ${OUT_DIR}/umi_deduplicated.bam: BAM file with deduplicated reads
# - ${OUT_DIR}/umi_tools_dedup.log: Log file for umi_tools deduplication

if [[ "$SKIP" -lt 5 ]]; then
    echo "Step 5a: Deduplicating 5' UMI BAM with umi_tools..."
    echo ""
    umi_tools dedup -I ${OUT_DIR}/umi_Aligned.sortedByCoord.out.bam \
                    --output-stats=${OUT_DIR}/umi_dedup_stats \
                    --paired \
                    --chimeric-pairs=discard \
                    --unpaired-reads=discard \
                    -S ${OUT_DIR}/umi_deduplicated.bam \
                    -L ${OUT_DIR}/umi_tools_dedup.log

    # --- INTERNAL READS ---
    # echo "Step 5b: Deduplicating Internal BAM with samtools..."
    # samtools collate -O -u -@ ${THREADS} ${OUT_DIR}/internal_Aligned.sortedByCoord.out.bam | \
    # samtools fixmate -m -u - - | \
    # samtools sort -u -@ ${THREADS} | \
    # samtools markdup -r -@ ${THREADS} - ${OUT_DIR}/internal_deduplicated.bam
    # ----------------------

    if [[ "$STEP" -eq 5 ]]; then echo "Stopping after Step 5."; exit 0; fi
else
    echo "Skipping Step 5..."
    echo ""
fi

# --- STEP 6: Final Production Indexing ---

# INPUTS:
# -I: Input BAM file with deduplicated reads

# OUTPUTS:
# - ${OUT_DIR}/umi_deduplicated.bam.bai: BAM index file

if [[ "$SKIP" -lt 6 ]]; then
    echo "Step 6: Indexing final deduplicated BAM files..."
    echo ""
    samtools index ${OUT_DIR}/umi_deduplicated.bam

    # --- INTERNAL READS ---
    # samtools index ${OUT_DIR}/internal_deduplicated.bam
    # ----------------------

    if [[ "$STEP" -eq 6 ]]; then echo "Stopping after Step 6."; exit 0; fi
else
    echo "Skipping Step 6..."
    echo ""
fi

# --- STEP 7: Gene Counting ---

# INPUTS:
# -p: Indicate that the input BAM contains paired-end reads
# -T: Number of threads to use for featureCounts
# -a: Path to GTF annotation file for gene features
# -t: Feature type to count (e.g., exon)
# -g: Attribute type to group features by (e.g., gene_id)
# -o: Output file for gene counts

if [[ "$SKIP" -lt 7 ]]; then
    echo "Step 7: Counting reads per gene with featureCounts..."
    echo ""
    featureCounts -p \
                  -T ${THREADS} \
                  -a d_data/refGenome/genomic.gtf \
                  -t exon \
                  -g gene_id \
                  -o ${OUT_DIR}/gene_counts.txt \
                  ${OUT_DIR}/umi_deduplicated.bam

    # --- INTERNAL READS ---
    # featureCounts -p \
    #               -T ${THREADS} \
    #               -a d_data/refGenome/genomic.gtf \
    #               -t exon \
    #               -g gene_id \
    #               -o ${OUT_DIR}/gene_counts.txt \
    #               ${OUT_DIR}/umi_deduplicated.bam ${OUT_DIR}/internal_deduplicated.bam
    # ----------------------


    if [[ "$STEP" -eq 7 ]]; then echo "Stopping after Step 7."; exit 0; fi
else
    echo "Skipping Step 7..."
    echo ""
fi

echo "========================================================"
echo "Pipeline complete! Final counts matrix generated."
echo "========================================================"