#!/bin/bash

cat SRR_Acc_List.txt | parallel -j 20 fastqc fastq/{}_1.fastq fastq/{}_2.fastq -o qc/raw
