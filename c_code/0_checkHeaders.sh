#!/bin/bash

OUTPUT_DIR="d_data/example_reads"
mkdir -p "$OUTPUT_DIR"

# Updated to use .fastq.gz based on your finding
for filepath in d_data/*.fastq.gz; do

    if [ -e "$filepath" ]; then
        filename=$(basename "$filepath")
        # Updated to strip .fastq.gz instead of _fastq.gz
        output_name="${filename%.fastq.gz}.txt"
        
        echo "========================================"
        echo "Processing: $filename"
        echo "========================================"
        
        # Switched from zcat to gzip -cd for better compatibility
        gzip -cd "$filepath" | head -n 40 | tee "$OUTPUT_DIR/$output_name"
    else
        echo "Warning: No matching FASTQ files found for '$filepath'"
        continue
    fi

done