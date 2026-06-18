#!/bin/bash

g++ -O3 -std=c++11 main.cpp -o cumi_tools -lz

# Define your massive input files here
R1="1_beads_1m_lib_S1_R1_001.fastq.gz"
R2="1_beads_1m_lib_S1_R2_001.fastq.gz"

echo "======================================================="
echo "Starting UMI Extraction Benchmark"
echo "Processing:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "======================================================="
echo "Running C++ Extractor..."
echo ""

# The 'time' command wraps the execution and measures the duration
time ./cumi_tools "$R1" "$R2"

echo ""
echo "======================================================="
echo "Job Complete!"
echo "======================================================="