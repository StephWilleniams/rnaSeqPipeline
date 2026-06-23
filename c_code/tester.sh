#!/bin/bash

# umi_tools extract \
#           --extract-method=regex \
#           --bc-pattern="(?P<discard_1>.{0,10})(?P<discard_2>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})(?P<discard_3>G{3,5})" \
#           -I d_data/3_beads_3m_lib_S3_R1_001.fastq.gz \
#           --read2-in=d_data/3_beads_3m_lib_S3_R2_001.fastq.gz \
#           -S o_outputs/extracted_R1.fastq.gz \
#           --read2-out=o_outputs/extracted_R2.fastq.gz \
#           --filtered-out=o_outputs/internal_R1.fastq.gz \
#           --filtered-out2=o_outputs/internal_R2.fastq.gz \
#           -L o_outputs/umi_tools_extract.log

./c_code/c_umi_extract -r d_data/3_beads_3m_lib_S3_R1_001.fastq.gz \
                    -R d_data/3_beads_3m_lib_S3_R2_001.fastq.gz \
                    -o o_outputs/extracted_R1.fastq.gz \
                    -O o_outputs/extracted_R2.fastq.gz \
                    -L o_outputs/c_umi_progress.log \
                    -i o_outputs/internal_R1.fastq.gz \
                    -I o_outputs/internal_R2.fastq.gz
