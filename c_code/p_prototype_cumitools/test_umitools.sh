#!/bin/bash

umi_tools extract \
          --extract-method=regex \
          --bc-pattern="(?P<discard_1>.*)(?P<discard_2>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})(?P<discard_3>G{3,5})" \
          -I test_R1.fastq.gz \
          --read2-in=test_R2.fastq.gz \
          -S umitest/extracted_R1.fastq.gz \
          --read2-out=umitest/extracted_R2.fastq.gz \
          --filtered-out=umitest/internal_R1.fastq.gz \
          --filtered-out2=umitest/internal_R2.fastq.gz \
          -L umitest/umi_tools_extract.log