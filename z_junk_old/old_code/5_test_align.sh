#!/bin/bash

STAR --runThreadN 8 \
     --genomeDir d_data/refGenome/mouse_star_index \
     --readFilesIn o_outputs/trimmed_fastqs/2_beads_2m_lib_S2_trimmed_R1.fastq o_outputs/trimmed_fastqs/2_beads_2m_lib_S2_trimmed_R2.fastq \
     --outSAMreadID Numbered \
     --outFileNamePrefix o_outputs/star_alignments/TEST_SINGLE_ \
     --outSAMtype BAM SortedByCoordinate