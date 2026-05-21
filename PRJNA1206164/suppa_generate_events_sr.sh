#!/bin/bash
###############################################################################
# suppa_generate_events_sr.sh
#
# Purpose : Generate alternative splicing events from the reference GTF using
#           SUPPA2, then compute per-event PSI values from Salmon TPM output.
#
# Inputs  : - Homo_sapiens.GRCh38.107_ERCC.gtf : Reference annotation
#           - transcript_tpm.txt                : Salmon TPM matrix
#                                                 (output of salmon_tx_tpm_counts.r)
#
# Outputs : - FC_events_<TYPE>_strict.ioe files : Splicing event definitions
#           - psi_<TYPE>.psi files               : Per-sample PSI values
#
# Next step: Run beta_regression_psi.r to test for differential splicing.
#
# Usage   : bash suppa_generate_events_sr.sh
###############################################################################

# ── Path configuration ────────────────────────────────────────────────────────
# ⚠️ Update these paths to match your local environment before running.
GTF_FILE="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"
TPM_FILE="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/transcript_tpm.txt"
OUTDIR="/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/suppa"

mkdir -p "${OUTDIR}"

echo "Starting SUPPA2 event generation at $(date)"

# ── Step 1: Generate splicing event definitions from GTF ─────────────────────
# Event types: SE (skipped exon), MX (mutually exclusive), SS (splice site),
#              RI (retained intron), FL (first/last exon)
suppa.py generateEvents \
    -i "${GTF_FILE}" \
    -o "${OUTDIR}/FC_events" \
    -f ioe \
    -e SE MX SS RI FL

echo "Event generation complete at $(date)"

# ── Step 2: Compute PSI per event for all 7 event types ──────────────────────
# Loops over each event type and produces a .psi file of per-sample PSI values
for event in SE RI MX AL AF A5 A3
do
    echo "  Computing PSI for event type: ${event}"
    suppa.py psiPerEvent \
        --ioe-file  "${OUTDIR}/FC_events_${event}_strict.ioe" \
        --expression-file "${TPM_FILE}" \
        -o "${OUTDIR}/psi_${event}"
done

echo "SUPPA2 PSI calculation complete at $(date)"
echo "PSI files written to: ${OUTDIR}"
