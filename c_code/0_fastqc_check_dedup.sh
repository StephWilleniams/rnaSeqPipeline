#!/bin/bash

# This script processes trimmed sequencing data by performing quality control using FastQC and MultiQC.

DIR1="o_outputs/qc_reports_dedup"
DIR2="o_outputs/qc_reports_dedup_multi"
if [ -d "$DIR1" ]; then
    echo "Directory '$DIR1' already exists. Deleting..."
    rm -rf "$DIR1"
fi
echo "Creating fresh directory '$DIR1'..."
mkdir -p "$DIR1"
if [ -d "$DIR2" ]; then
    echo "Directory '$DIR2' already exists. Deleting..."
    rm -rf "$DIR2"
fi
echo "Creating fresh directory '$DIR2'..."
mkdir -p "$DIR2"

for file in o_outputs/sample_*/umi_deduplicated.bam; do
    # Extract the sample name (e.g., "sample_A") from the path
    sample_name=$(basename "$(dirname "$file")")
    
    echo "Processing $sample_name..."
    
    # Run FastQC on the single file
    fastqc -o o_outputs/qc_reports_dedup/ -t 1 "$file"
    
    # Rename the resulting outputs to include the sample name
    mv o_outputs/qc_reports_dedup/umi_deduplicated_fastqc.html o_outputs/qc_reports_dedup/${sample_name}_dedup.html
    mv o_outputs/qc_reports_dedup/umi_deduplicated_fastqc.zip o_outputs/qc_reports_dedup/${sample_name}_dedup.zip
done

# ==============================================================================
# MultiQC v1.35 Workaround for BAM files on Google Drive CloudStorage
# ==============================================================================

echo "Setting up local scratch directory for MultiQC processing..."
TMP_QC_DIR=$(mktemp -d /tmp/multiqc_scratch.XXXXXX)

echo "Extracting raw fastqc_data.txt logs to bypass filesystem and zip stream limits..."
for zipfile in o_outputs/qc_reports_dedup/*_dedup.zip; do
    if [ -f "$zipfile" ]; then
        bname=$(basename "$zipfile" .zip)
        unzip -p "$zipfile" "*/fastqc_data.txt" > "$TMP_QC_DIR/${bname}_fastqc_data.txt" 2>/dev/null
    fi
done

echo "Running MultiQC locally on extracted text files..."
# Run multiqc inside the temp directory so it outputs natively there
(
    cd "$TMP_QC_DIR" || exit 1
    multiqc --fn_as_s_name --force .
)

echo "Moving generated MultiQC reports back to $DIR2..."
if [ -f "$TMP_QC_DIR/multiqc_report.html" ]; then
    mv "$TMP_QC_DIR/multiqc_report.html" "$DIR2/"
    if [ -d "$TMP_QC_DIR/multiqc_data" ]; then
        mv "$TMP_QC_DIR/multiqc_data" "$DIR2/"
    fi
    echo "MultiQC processing completed successfully."
else
    echo "Error: MultiQC failed to generate a report."
fi

echo "Cleaning up local scratch directory..."
rm -rf "$TMP_QC_DIR"