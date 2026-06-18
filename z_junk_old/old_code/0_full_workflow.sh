#!/bin/bash
set -e

# --- CONFIGURATION ---
THREADS=8
R1_IN="d_data/1_beads_1m_lib_S1_R1_001.fastq.gz"
R2_IN="d_data/1_beads_1m_lib_S1_R2_001.fastq.gz"
STAR_INDEX="d_data/refGenome/mouse_star_index/"
OUT_DIR="processed_data"
mkdir -p ${OUT_DIR}

echo "========================================================"
echo "Starting Verified Takara SMART-Seq Pipeline (STAR)"
echo "========================================================"
echo ""

# --- STEP 1: fastp (Quality, Adapter Cleanup, and minimum string length) ---
echo "Step 1: Running fastp..."

fastp -i ${R1_IN} -I ${R2_IN} \
      -o ${OUT_DIR}/trimmed_R1.fastq.gz -O ${OUT_DIR}/trimmed_R2.fastq.gz \
      --html ${OUT_DIR}/fastp_report.html --json ${OUT_DIR}/fastp_report.json \
      -q 20 -u 30 \
      --thread ${THREADS} 2>${OUT_DIR}/fastp_run.log # Redirecting stderr to a log file for fastp

# --- STEP 2: Official UMI Extraction (Corrected Anchor Pattern) ---
echo "Step 2: Extracting UMIs..."

# (?P<discard_3>G{3,5})

umi_tools extract \
      --extract-method=regex \
      --bc-pattern="(?P<discard_1>.)(?P<discard_2>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})" \
      -I ${OUT_DIR}/trimmed_R1.fastq.gz \
      --read2-in=${OUT_DIR}/trimmed_R2.fastq.gz \
      -S ${OUT_DIR}/extracted_R1.fastq.gz \
      --read2-out=${OUT_DIR}/extracted_R2.fastq.gz \
      --filtered-out=${OUT_DIR}/internal_R1.fastq.gz \
      --filtered-out2=${OUT_DIR}/internal_R2.fastq.gz \
      -L ${OUT_DIR}/umi_tools_extract.log

# --- STEP 3: Alignment (STAR Native Compressed Multi-File Stream) ---
echo "Step 3: Aligning compiled compressed fractions with STAR..."

# CONDA_SUBDIR=osx-64 mamba install -y bioconda::star=2.7.10b

STAR --runThreadN ${THREADS} \
     --genomeDir ${STAR_INDEX} \
     --readFilesIn ${OUT_DIR}/extracted_R1.fastq.gz,${OUT_DIR}/internal_R1.fastq.gz \
                   ${OUT_DIR}/extracted_R2.fastq.gz,${OUT_DIR}/internal_R2.fastq.gz \
     --readFilesCommand "gzip -d -c" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes Standard \
     --outFileNamePrefix ${OUT_DIR}/star_

# conda clean --all -y
# conda remove --force star -y

# --- STEP 4: Indexing (samtools) ---
echo "Step 4: Indexing BAM file..."

samtools index ${OUT_DIR}/star_Aligned.sortedByCoord.out.bam

# --- STEP 5: Deduplication ---
echo "Step 5: Deduplicating BAM file with UMICollapse..."

umi_tools dedup -I ${OUT_DIR}/star_Aligned.sortedByCoord.out.bam \
                --output-stats=${OUT_DIR}/dedup_stats \
                --paired \
                -S ${OUT_DIR}/final_deduplicated.bam

# --- STEP 6: Final Production Indexing ---
echo "Step 6: Indexing final deduplicated BAM file..."

samtools index ${OUT_DIR}/final_deduplicated.bam

# --- STEP 7: Gene Counting (featureCounts) ---
echo "Step 7: Counting reads per gene with featureCounts..."

featureCounts -p -T ${THREADS} \
              -a d_data/refGenome/genomic.gtf \
              -t exon \
              -g gene_id \
              -o ${OUT_DIR}/gene_counts.txt \
              ${OUT_DIR}/final_deduplicated.bam

echo "========================================================"
echo "Pipeline complete! Final file: ${OUT_DIR}/final_deduplicated.bam"
echo "========================================================"

