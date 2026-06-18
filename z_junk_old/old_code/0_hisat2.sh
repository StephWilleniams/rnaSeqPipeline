#!/bin/bash

# Extract splice sites from the GTF
python /opt/miniconda3/envs/umi_mapping/bin/extract_splice_sites.py d_data/refGenome/genomic.gtf > d_data/refGenome/mouse.ss

# Extract exons from the GTF
python /opt/miniconda3/envs/umi_mapping/bin/extract_exons.py d_data/refGenome/genomic.gtf > d_data/refGenome/mouse.exon

hisat2-build -p 2 \
             --ss d_data/refGenome/mouse.ss \
             --exon d_data/refGenome/mouse.exon \
             d_data/refGenome/genomic.fna \
             d_data/refGenome/mouse_hisat2_index
