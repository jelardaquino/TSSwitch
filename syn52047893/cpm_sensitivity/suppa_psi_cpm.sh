#!/usr/bin/env bash
# PSI per event from CPM instead of TPM. Same IOE event definitions as the
# TPM-based run (suppa_bambu_gtf/), so only the expression input changes.
set -euo pipefail

BASE=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis
SUPPA_ENV=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/jaquino_conda_env/r4_5
OUT_DIR=$BASE/suppa_bambu_gtf_cpm

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

Rscript "$BASE/scripts/tx_counts_to_cpm.r"

for event in SE RI MX AL AF A5 A3; do
    echo "psiPerEvent: $event"
    "$SUPPA_ENV/bin/python" "$SUPPA_ENV/bin/suppa.py" psiPerEvent \
        --ioe-file "$BASE/suppa_bambu_gtf/FC_events_${event}_strict.ioe" \
        --expression-file "$BASE/bambu_tx_cpm.txt" \
        -o "$OUT_DIR/psi_${event}"
done

echo "Done: $OUT_DIR"
