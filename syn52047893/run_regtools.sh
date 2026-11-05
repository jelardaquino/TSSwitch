#!/bin/bash

cat samples.txt | parallel -j 12 'regtools junctions extract -a 8 -m 50 -M 500000 -s XS -o {}_junctions.bed ../{}.bam'
