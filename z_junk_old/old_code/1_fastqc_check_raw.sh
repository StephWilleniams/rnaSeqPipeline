#!/bin/bash

# This script processes raw sequencing data by performing quality control using FastQC and MultiQC.

DIR="o_outputs/qc_reports"
if [ -d "$DIR" ]; then
    echo "Directory '$DIR' already exists. Deleting..."
    rm -rf "$DIR"
fi
echo "Creating fresh directory '$DIR'..."
mkdir -p "$DIR"

fastqc -o o_outputs/qc_reports/ -t 12 d_data/*.fastq.gz

multiqc -o o_outputs/qc_reports_multi/ o_outputs/qc_reports/
