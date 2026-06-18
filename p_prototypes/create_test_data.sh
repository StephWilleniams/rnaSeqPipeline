#!/bin/bash

# Convert a subset of 1_beads_1m_lib_S1_R1_001.fastq.gz into test_R1.fastq.gz
# and the corresponding subset of 1_beads_1m_lib_S1_R2_001.fastq.gz into test_R2.fastq.gz
gzcat 1_beads_1m_lib_S1_R1_001.fastq.gz | head -n 4000 | gzip > test_R1.fastq.gz
gzcat 1_beads_1m_lib_S1_R2_001.fastq.gz | head -n 4000 | gzip > test_R2.fastq.gz
