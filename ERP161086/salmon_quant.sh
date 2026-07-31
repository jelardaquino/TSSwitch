#!/bin/bash
echo "Starting Salmon quant at $(date)"
cat /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/SRR_Acc_List.txt | parallel -j 10 salmon quant -i /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn21589959/scripts/hs_cdna_index -l A -1 /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/fastq/{}_1.fastq.gz -2 /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/fastq/{}_2.fastq.gz -o /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/salmon/{}_quant --gcBias --validateMappings -p 8
echo "Finished Salmon quant at $(date)"
