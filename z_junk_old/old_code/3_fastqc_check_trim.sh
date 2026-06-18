#!/bin/bash

# This script processes trimmed sequencing data by performing quality control using FastQC and MultiQC.

DIR1="o_outputs/qc_reports_trimmed"
DIR2="o_outputs/qc_reports_trimmed_multi"
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

fastqc -o o_outputs/qc_reports_trimmed/ -t 12 o_outputs/trimmed_fastqs/*.fastq.gz

multiqc -o o_outputs/qc_reports_trimmed_multi/ o_outputs/qc_reports_trimmed/