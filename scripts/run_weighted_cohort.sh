#!/usr/bin/env bash
# Score, weighted prepare and infer for every sample in a directory, plus an
# unweighted baseline for comparison.
#
# Layout expected:   <samples_dir>/<sample>/aligned_to_13v2.bam   BAM aligned to chm13v2 (input to prepare)
#                    <samples_dir>/<sample>/extracted.tsv         modkit extract read-calls table on hg38 (input to score)
# Output:            <out_dir>/<sample>/scores/    read_scores.parquet, moments.json, score_summary.json
#                    <out_dir>/<sample>/weighted/  methylation.bed + predictions with read weights
#                    <out_dir>/<sample>/baseline/  methylation.bed + predictions without weights
#                    <out_dir>/summary.csv         top prediction per sample, baseline vs weighted
#
# Usage:
#   scripts/run_weighted_cohort.sh <samples_dir> <atlas.feather> <out_dir> [threads]
#
# Environment overrides:
#   BAM_NAME       BAM file name inside each sample dir      (default aligned_to_13v2.bam)
#   CALLS_NAME     read-calls table inside each sample dir   (default extracted.tsv)
#   SCORE_ARGS     extra options for `nanoflux score`, e.g. "--weight-temperature 0.5 --weight-min 0.1"
#   PREPARE_ARGS   extra options for `nanoflux prepare`      (default "--skip-alignment")
#   RUN_BASELINE   1 to also run the unweighted pipeline     (default 1)
#   CLEANUP        1 to delete the sorted BAM copies prepare writes into the output dirs (default 1)
#   SAMPLES        space-separated subset of sample names    (default: every subdirectory)
#
# Needs the nanoflux environment active (nanoflux, samtools, modkit, bedtools on PATH),
# installed from this branch: `uv pip install -e .` or `pip install -e .`.
# Steps whose outputs already exist are skipped, so the script can be re-run.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '2,25p' "$0"
    exit 1
fi

SAMPLES_DIR=$(cd "$1" && pwd)
ATLAS=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
OUT=$3
THREADS=${4:-8}
BAM_NAME=${BAM_NAME:-aligned_to_13v2.bam}
CALLS_NAME=${CALLS_NAME:-extracted.tsv}
SCORE_ARGS=${SCORE_ARGS:-}
PREPARE_ARGS=${PREPARE_ARGS:-"--skip-alignment"}
RUN_BASELINE=${RUN_BASELINE:-1}
CLEANUP=${CLEANUP:-1}

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
LOG="$OUT/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

for tool in nanoflux samtools modkit bedtools; do
    command -v "$tool" >/dev/null || { echo "[ERROR] $tool not on PATH; activate the nanoflux environment" >&2; exit 1; }
done
[[ -f "$ATLAS" ]] || { echo "[ERROR] atlas not found: $ATLAS" >&2; exit 1; }

if [[ -n "${SAMPLES:-}" ]]; then
    read -r -a sample_list <<< "$SAMPLES"
else
    sample_list=()
    for d in "$SAMPLES_DIR"/*/; do sample_list+=("$(basename "$d")"); done
fi
echo "[INFO] $(date)  ${#sample_list[@]} samples, atlas $ATLAS, out $OUT, threads $THREADS"
echo "[INFO] score args: '${SCORE_ARGS}'  prepare args: '${PREPARE_ARGS}'  baseline: $RUN_BASELINE"

run_pipeline() {  # <sample> <bam> <run_dir> [extra prepare args...]
    local s=$1 bam=$2 run_dir=$3
    shift 3
    if [[ -f "$run_dir/predicted_labels.json" ]]; then
        echo "[SKIP]  $s $(basename "$run_dir"): predictions exist"
        return 0
    fi
    rm -rf "$run_dir"
    echo "[START] $s $(basename "$run_dir") prepare"
    # shellcheck disable=SC2086
    nanoflux prepare -i "$bam" -o "$run_dir" -c -t "$THREADS" $PREPARE_ARGS "$@"
    echo "[START] $s $(basename "$run_dir") infer"
    nanoflux infer -i "$run_dir/methylation.bed" -o "$run_dir" -c -t "$THREADS"
    if [[ "$CLEANUP" == 1 ]]; then
        rm -f "$run_dir"/aligned_to_*.bam "$run_dir"/aligned_to_*.bam.bai "$run_dir"/aligned_to_*.sam
    fi
    echo "[DONE]  $s $(basename "$run_dir")"
}

failed=()
for s in "${sample_list[@]}"; do
    dir="$SAMPLES_DIR/$s"
    bam="$dir/$BAM_NAME"
    calls="$dir/$CALLS_NAME"
    if [[ ! -f "$bam" || ! -f "$calls" ]]; then
        echo "[SKIP]  $s: missing $BAM_NAME or $CALLS_NAME"
        continue
    fi
    o="$OUT/$s"
    (
        set -e
        # 1. score reads against the atlas (hg38 read table)
        if [[ -f "$o/scores/read_scores.parquet" ]]; then
            echo "[SKIP]  $s score: exists"
        else
            rm -rf "$o/scores"
            echo "[START] $s score"
            # shellcheck disable=SC2086
            nanoflux score -i "$calls" --atlas "$ATLAS" -o "$o/scores" -c $SCORE_ARGS
        fi
        # 2. weighted prepare + infer (chm13v2 BAM, weights joined by read name)
        run_pipeline "$s" "$bam" "$o/weighted" --read-weights "$o/scores/read_scores.parquet"
        # 3. unweighted baseline
        if [[ "$RUN_BASELINE" == 1 ]]; then
            run_pipeline "$s" "$bam" "$o/baseline"
        fi
    ) || { echo "[FAILED] $s"; failed+=("$s"); }
done

# 4. summary table: top prediction per sample and run
python - "$OUT" <<'EOF'
import json, sys
from pathlib import Path
import pandas as pd

out = Path(sys.argv[1])
rows = []
for sample_dir in sorted(p for p in out.iterdir() if p.is_dir()):
    for run in ("baseline", "weighted"):
        probs = sample_dir / run / "predicted_probabilities.csv"
        if not probs.exists():
            continue
        p = pd.read_csv(probs, index_col=0)
        top = p.iloc[0]
        row = {"sample": sample_dir.name, "run": run, "top_class": p.index[0], "top_prob": float(top.iloc[0])}
        info = sample_dir / run / "weighted_pileup_info.json"
        if info.exists():
            d = json.loads(info.read_text())
            row.update({"probes_covered": d.get("n_probes_covered"), "reads_unscored": d.get("n_reads_unscored")})
        rows.append(row)
if rows:
    df = pd.DataFrame(rows).sort_values(["sample", "run"])
    df.to_csv(out / "summary.csv", index=False)
    print(df.to_string(index=False))
EOF

echo "[INFO] $(date)  finished; log: $LOG"
if [[ ${#failed[@]} -gt 0 ]]; then
    echo "[INFO] failed samples: ${failed[*]}"
    exit 1
fi
