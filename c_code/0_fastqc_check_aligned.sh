#!/bin/bash

# This script processes trimmed sequencing data by performing quality control using FastQC and MultiQC.

DIR1="o_outputs/qc_reports_aligned"
DIR2="o_outputs/qc_reports_aligned_multi"
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

for file in o_outputs/sample_*/umi_Aligned.sortedByCoord.out.bam; do
    # Extract the sample name (e.g., "sample_A") from the path
    sample_name=$(basename "$(dirname "$file")")
    
    echo "Processing $sample_name..."
    
    # Run FastQC on the single file
    fastqc -o o_outputs/qc_reports_aligned/ -t 1 "$file"
    
    # Rename the resulting outputs to include the sample name
    mv o_outputs/qc_reports_aligned/umi_Aligned.sortedByCoord.out_fastqc.html o_outputs/qc_reports_aligned/${sample_name}_aligned.html
    mv o_outputs/qc_reports_aligned/umi_Aligned.sortedByCoord.out_fastqc.zip o_outputs/qc_reports_aligned/${sample_name}_aligned.zip
done

multiqc --fn_as_s_name --force -o o_outputs/qc_reports_aligned_multi/ o_outputs/qc_reports_aligned/