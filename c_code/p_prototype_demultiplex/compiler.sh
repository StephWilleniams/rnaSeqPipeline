#!/bin/bash

g++ -O3 -std=c++11 main.cpp -o c_custom_demultiplex -lz

./c_code/p_prototype_demultiplex/c_custom_demultiplex -r d_data/Undetermined_S0_R1_001.fastq.gz -R d_data/Undetermined_S0_R2_001.fastq.gz -o d_data/bonusReads/ -m 3
