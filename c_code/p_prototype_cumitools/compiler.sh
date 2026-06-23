#!/bin/bash

# g++ -O3 -std=c++11 main.cpp edlib.cpp -o cumi_tools -lz # edlib

g++ -O3 -std=c++11 main.cpp -o c_umi_extract -lz

# ./c_umi_extract -r 1_beads_1m_lib_S1_R1_001.fastq.gz \
#              -R 1_beads_1m_lib_S1_R2_001.fastq.gz \
#              -o extracted_R1.fastq.gz \
#              -O extracted_R2.fastq.gz \
#              -L progress.log \
#              -i internal_R1.fastq.gz \
#              -I internal_R2.fastq.gz
