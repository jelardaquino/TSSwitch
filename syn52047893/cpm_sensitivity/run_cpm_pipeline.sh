#!/usr/bin/env bash
# Runs the CPM pipeline after step 1 (DIU). Steps 2a, 2b and the promoter chain
# (03 -> 04 -> 05) run in parallel; step 06 compares everything with the TPM results.
# Usage: bash run_cpm_pipeline.sh [--with-diu]
set -uo pipefail

J=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis
S=$J/scripts/cpm_pipeline
C=$J/cpm_pipeline
LOG=$C/logs

run_step() {
    name=$1; script=$2; workdir=$3
    echo "[$(date '+%F %T')] start $name"
    (cd "$workdir" && Rscript "$script" > "$LOG/$name.log" 2>&1)
    status=$?
    echo "[$(date '+%F %T')] end   $name (exit $status)"
    return $status
}

if [ "${1:-}" = "--with-diu" ]; then
    run_step 01_diu "$S/01_diu_dexseq_cpm.r" "$C/diu" || exit 1
fi

[ -s "$C/diu/sig_transcripts.txt" ] || { echo "Step 1 outputs missing"; exit 1; }

run_step 02a_psi_additive "$S/02a_psi_betareg_additive_cpm.r" "$C/psi_additive" &
pid_a=$!
run_step 02b_psi_condition_x_sex "$S/02b_psi_betareg_condition_x_sex_cpm.r" "$C/psi_condition_x_sex" &
pid_b=$!
(
    run_step 03_promoter_usage "$S/03_promoter_usage_cpm.r" "$C/promoter_usage" &&
    run_step 04_motif_tf "$S/04_motif_tf_cpm.r" "$C/promoter_usage" &&
    run_step 05_promoter_isoform_merge "$S/05_promoter_isoform_merge_cpm.r" "$C/promoter_usage"
) &
pid_p=$!

fail=0
for p in $pid_a $pid_b $pid_p; do wait $p || fail=1; done

run_step 06_compare "$S/06_compare_tpm_vs_cpm.r" "$C/comparison" || fail=1
exit $fail
