#!/bin/bash
echo "Starting Salmon quant at $(date)"
cat SRR_Acc_List.txt | parallel -j 10 salmon quant -i /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn21589959/scripts/hs_cdna_index -l A -1 trimmed/{}_R1_paired.fastq -2 trimmed/{}_R2_paired.fastq -o {}_quant --gcBias --validateMappings -p 8
echo "Finished Salmon quant at $(date)"
