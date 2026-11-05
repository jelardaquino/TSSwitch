#!/usr/bin/env bash
# Quantify the ERP161086 short-read samples against the long-read Bambu extended
# annotation (known + novel transcripts). Index and quant settings mirror the
# Ensembl 107 run (Reference/salmon_indexv107/MakeSalmonIndex.sh and
# ERP161086_quantSummary/salmon_quantv2): gffread transcript FASTA, k = 31,
# no decoys, Salmon 1.10.3, -l A --gcBias --validateMappings.
set -euo pipefail

S=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing
ENV=$S/jaquino_conda_env
GFFREAD=$ENV/rnaseq/bin/gffread
SALMON=$ENV/salmon_1103/bin/salmon

GTF=$S/syn52065646_LongRead/jaquino_analysis/NanoporeRNASeq_extended_annotations.gtf
GENOME=$S/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38_ERCC.fa
FASTQ_DIR=$S/ERP161086/fastq
METADATA=$S/ERP161086/metadata.txt

OUT=$S/ERP161086/salmon_bambu
INDEX=$OUT/index/salmon_index_bambu_extended
QUANT=$OUT/quant

INDEX_THREADS=32
PARALLEL_SAMPLES=4
THREADS_PER_SAMPLE=16

mkdir -p "$OUT/index" "$QUANT" "$OUT/logs"
cd "$OUT"

if [ ! -s bambu_extended_transcripts.fa ]; then
    echo "[$(date)] gffread: building transcript FASTA"
    "$GFFREAD" "$GTF" -g "$GENOME" -w bambu_extended_transcripts.fa
fi
echo "Transcripts in FASTA: $(grep -c '^>' bambu_extended_transcripts.fa)"

if [ ! -f "$INDEX/versionInfo.json" ]; then
    echo "[$(date)] salmon index"
    "$SALMON" index -t bambu_extended_transcripts.fa -i "$INDEX" -k 31 -p "$INDEX_THREADS" \
        > "$OUT/logs/salmon_index.log" 2>&1
fi

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

tail -n +2 "$METADATA" | cut -f1 | xargs -P "$PARALLEL_SAMPLES" -I{} bash -c 'quant_one {}'

echo "[$(date)] All samples quantified: $(ls -d "$QUANT"/*_quant | wc -l)"
