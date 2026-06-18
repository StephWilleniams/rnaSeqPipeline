#!/bin/bash

# This script generates a STAR index for the mouse reference genome, which is necessary for aligning RNA-seq reads to the genome.

DIR="d_data/refGenome/mouse_star_index"
if [ -d "$DIR" ]; then
    echo "Directory '$DIR' already exists. Deleting..."
    rm -rf "$DIR"
fi
echo "Creating fresh directory '$DIR'..."
mkdir -p "$DIR"

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir d_data/refGenome/mouse_star_index \
     --genomeFastaFiles d_data/refGenome/GCF_000001635.27_GRCm39_genomic.fna \
     --sjdbGTFfile d_data/refGenome/GCF_000001635.27_GRCm39_genomic.gff \
     --sjdbOverhang 99
