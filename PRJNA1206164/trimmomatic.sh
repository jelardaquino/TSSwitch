#!/bin/bash

cat SRR_Acc_List.txt | \
parallel -j 10 '
  trimmomatic PE \
    -threads 8 \
    fastq/{}_1.fastq fastq/{}_2.fastq \
    trimmed/{}_R1_paired.fastq trimmed/{}_R1_unpaired.fastq \
    trimmed/{}_R2_paired.fastq trimmed/{}_R2_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36
'
