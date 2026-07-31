#!/bin/bash

cat /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/SRR_Acc_List.txt | parallel -j 14 fastqc /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/fastq/{}_1.fastq.gz /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/fastq/{}_2.fastq.gz -o /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/qc/raw
