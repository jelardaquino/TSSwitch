#!/usr/bin/env bash
# Quantify all 41 E-MTAB-14167 prefrontal cortex samples against the long-read
# Bambu extended annotation, reusing the index built for the ERP161086 7 vs 7 run
# (ERP161086/scripts/salmon_bambu_index_quant.sh). Same settings: Salmon 1.10.3,
# raw FASTQ, -l A --gcBias --validateMappings.
set -euo pipefail

S=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing
SALMON=$S/jaquino_conda_env/salmon_1103/bin/salmon
INDEX=$S/ERP161086/salmon_bambu/index/salmon_index_bambu_extended
FASTQ_DIR=$S/E-MTAB-14167/Fastq
METADATA=$S/E-MTAB-14167/Step01_metadata/sample_metadata_matchedv2.tsv

OUT=$S/ERP161086/emtab14167_replication/salmon_bambu
QUANT=$OUT/quant

PARALLEL_SAMPLES=6
THREADS_PER_SAMPLE=12

mkdir -p "$QUANT" "$OUT/logs"

quant_one() {
    run=$1
    out=$QUANT/${run}_quant
    if [ -s "$out/quant.sf" ]; then
        echo "[$(date)] $run already quantified, skipping"
        return 0
    fi
    echo "[$(date)] salmon quant $run"
    "$SALMON" quant -i "$INDEX" -l A \
        -1 "$FASTQ_DIR/${run}_1.fastq.gz" -2 "$FASTQ_DIR/${run}_2.fastq.gz" \
        --gcBias --validateMappings -p "$THREADS_PER_SAMPLE" \
        -o "$out" > "$OUT/logs/salmon_quant_${run}.log" 2>&1
    echo "[$(date)] done $run: mapping rate $(grep -o '"percent_mapped": [0-9.]*' "$out/aux_info/meta_info.json")"
}
export -f quant_one
export SALMON INDEX QUANT FASTQ_DIR THREADS_PER_SAMPLE OUT

# PFC runs only (column 8 = ENA_RUN, column 11 = RegionShort)
awk -F'\t' 'NR > 1 && $11 == "PFC" {print $8}' "$METADATA" | \
    xargs -P "$PARALLEL_SAMPLES" -I{} bash -c 'quant_one {}'

echo "[$(date)] All samples quantified: $(ls -d "$QUANT"/*_quant | wc -l)"
