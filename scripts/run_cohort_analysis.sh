#!/usr/bin/env bash
# Whole analysis for a labelled cohort: prepare, score against the control atlas, build tumour
# models (leave-one-out for the samples they are built from), score against them, and compare
# all approaches in one table.
#
# Usage:
#   scripts/run_cohort_analysis.sh <sheet.tsv> <control_atlas.feather> <work_dir> [threads]
#
# Sheet: tab-separated with a header line and the columns
#   sample_id    name, used for the output folders
#   category     control | tumor_low | tumor_high     (tumor_high samples build the tumour models)
#   true_label   class exactly as the classifier names it, e.g. "CONTR, INFLAM", "GBM", "MB, G3 G4"
#   bam          BAM aligned to chm13v2 with MM/ML tags            (input of nanoflux prepare)
#   calls        modkit read-calls table of the sample on hg38     (input of nanoflux score)
#
# Steps (each is skipped when its output exists, so the script can be re-run):
#   1  <work>/<sample>/baseline/       nanoflux prepare --extract-reads   -> reads.tsv, methylation.bed
#   2  <work>/<sample>/scores/         nanoflux score                     -> z, llr_per_call
#   3  <work>/tumor_atlas/full/        tumour atlas (+ Beta prior, + calibration) from all tumor_high samples
#      <work>/tumor_atlas/without_<s>/ the same without sample s, for every tumor_high sample (unless BIASED=1)
#   4  <work>/<sample>/scores_tumor/   nanoflux score --tumor-atlas ...   -> all five scores: z_control,
#                                       llr_control_random_per_call, z_tumor, llr_tumor_random_per_call,
#                                       llr_tumor_control_per_call
#   5  <work>/results/                 scripts/evaluate_approaches.py     -> approaches_table.txt,
#                                       approaches_summary.csv, approaches_long.csv, tumor_low_*.csv
#
# Environment overrides:
#   PREPARE_ARGS            extra options for nanoflux prepare    (default "--skip-alignment")
#   SCORE_ARGS              extra options for nanoflux score      (default none)
#   MOMENTS_FROM_CONTROLS   1: refit the caller moments on the control samples of the sheet instead of
#                           using the ones shipped with the package                       (default 0)
#   TUMOR_MODEL             prior: Beta prior on the tumour atlas (--tumor-prior auto, deep atlas);
#                           calibration: calibration-table model (sparse atlas)           (default prior)
#   TUMOR_PRIOR_MIN_COV     coverage for the tumour prior fit, build and score            (default 20)
#   BIASED                  1: no leave-one-out; every sample, including the pooled tumor_high ones, is
#                           scored against tumor_atlas/full (in-pool, optimistic for tumor_high) (default 0)
#   CLEANUP                 1: delete the sorted BAM copies prepare writes                (default 1)
#   PYTHON                  python of the environment nanoflux is installed in            (default python)
#
# Needs nanoflux (this branch, `pip install -e .`), samtools, modkit 0.2.5 and bedtools on PATH.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '2,33p' "$0"
    exit 1
fi

SHEET=$1
ATLAS=$2
WORK=$3
THREADS=${4:-8}
PREPARE_ARGS=${PREPARE_ARGS:-"--skip-alignment"}
SCORE_ARGS=${SCORE_ARGS:-}
MOMENTS_FROM_CONTROLS=${MOMENTS_FROM_CONTROLS:-0}
TUMOR_MODEL=${TUMOR_MODEL:-prior}
TUMOR_PRIOR_MIN_COV=${TUMOR_PRIOR_MIN_COV:-20}
BIASED=${BIASED:-0}
CLEANUP=${CLEANUP:-1}
PYTHON=${PYTHON:-python}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for tool in nanoflux samtools modkit bedtools "$PYTHON"; do
    command -v "$tool" >/dev/null || { echo "[ERROR] $tool not on PATH" >&2; exit 1; }
done
[[ -f "$SHEET" ]] || { echo "[ERROR] sheet not found: $SHEET" >&2; exit 1; }
[[ -f "$ATLAS" ]] || { echo "[ERROR] atlas not found: $ATLAS" >&2; exit 1; }
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
exec > >(tee -a "$WORK/run_$(date +%Y%m%d_%H%M%S).log") 2>&1
echo "[INFO] $(date)  sheet $SHEET  atlas $ATLAS  work $WORK  threads $THREADS"

# ---- read the sheet (columns located by name) ----
header=$(head -1 "$SHEET")
col() { echo "$header" | tr '\t' '\n' | grep -n -x "$1" | cut -d: -f1; }
for c in sample_id category true_label bam calls; do
    [[ -n "$(col $c)" ]] || { echo "[ERROR] sheet has no column '$c'" >&2; exit 1; }
done
field() { awk -F'\t' -v r="$1" -v c="$(col "$2")" 'NR == r + 1 { print $c }' "$SHEET"; }
n_samples=$(( $(grep -c . "$SHEET") - 1 ))

samples=(); categories=(); bams=(); calls=()
for i in $(seq 1 "$n_samples"); do
    samples+=("$(field "$i" sample_id)"); categories+=("$(field "$i" category)")
    bams+=("$(field "$i" bam)"); calls+=("$(field "$i" calls)")
    [[ -f "${bams[$((i-1))]}" ]] || { echo "[ERROR] BAM missing: ${bams[$((i-1))]}" >&2; exit 1; }
    [[ -f "${calls[$((i-1))]}" ]] || { echo "[ERROR] calls missing: ${calls[$((i-1))]}" >&2; exit 1; }
done
echo "[INFO] $n_samples samples: ${samples[*]}"

# ---- optional: caller moments from the controls of this cohort ----
moments_arg=()
if [[ "$MOMENTS_FROM_CONTROLS" == 1 ]]; then
    if [[ ! -f "$WORK/moments/moments.json" ]]; then
        controls=()
        for i in "${!samples[@]}"; do [[ "${categories[$i]}" == control ]] && controls+=("${calls[$i]}"); done
        echo "[START] caller moments from ${#controls[@]} controls"
        rm -rf "$WORK/moments"
        nanoflux score -i "${controls[0]}" --atlas "$ATLAS" -o "$WORK/moments" -c --moments-from "${controls[@]}"
    fi
    moments_arg=(--moments "$WORK/moments/moments.json")
fi

# ---- 1 + 2: prepare and score per sample ----
for i in "${!samples[@]}"; do
    s=${samples[$i]}; o="$WORK/$s"; mkdir -p "$o"
    if [[ -f "$o/baseline/reads.tsv" && -f "$o/baseline/methylation.bed" ]]; then
        echo "[SKIP]  $s prepare"
    else
        rm -rf "$o/baseline"
        echo "[START] $s prepare"
        # shellcheck disable=SC2086
        nanoflux prepare -i "${bams[$i]}" -o "$o/baseline" -c -t "$THREADS" $PREPARE_ARGS --extract-reads
        [[ "$CLEANUP" == 1 ]] && rm -f "$o"/baseline/aligned_to_*.bam "$o"/baseline/aligned_to_*.bam.bai "$o"/baseline/aligned_to_*.sam
    fi
    if [[ -f "$o/scores/read_scores.parquet" ]]; then
        echo "[SKIP]  $s score"
    else
        rm -rf "$o/scores"
        echo "[START] $s score"
        # shellcheck disable=SC2086
        nanoflux score -i "${calls[$i]}" --atlas "$ATLAS" -o "$o/scores" -c ${moments_arg[@]+"${moments_arg[@]}"} $SCORE_ARGS
    fi
done

# ---- 3: tumour models ----
high_idx=()
for i in "${!samples[@]}"; do [[ "${categories[$i]}" == tumor_high ]] && high_idx+=("$i"); done
build_model() {  # <name> <index to leave out, or -1>
    local name=$1 leave=$2 out="$WORK/tumor_atlas/$1" files=()
    if [[ -f "$out/build_info.json" ]]; then echo "[SKIP]  tumour model $name"; return 0; fi
    for i in "${high_idx[@]}"; do [[ "$i" != "$leave" ]] && files+=("${calls[$i]}"); done
    echo "[START] tumour model $name from ${#files[@]} samples"
    rm -rf "$out"
    "$PYTHON" "$HERE/build_tumor_atlas.py" --control-atlas "$ATLAS" --samples "${files[@]}" --out "$out" \
        --prior-min-cov "$TUMOR_PRIOR_MIN_COV" --allow-shallow-prior
}
if [[ ${#high_idx[@]} -ge 2 ]]; then
    build_model full -1
    if [[ "$BIASED" == 1 ]]; then
        echo "[INFO]  BIASED=1: no leave-one-out models, pooled tumor_high samples are scored against full/"
    else
        for i in "${high_idx[@]}"; do build_model "without_${samples[$i]}" "$i"; done
    fi
    fit_ok=$("$PYTHON" -c "import json,sys; print(json.load(open(sys.argv[1])).get('prior',{}).get('fit_ok','unknown'))" \
        "$WORK/tumor_atlas/full/build_info.json" 2>/dev/null || echo unknown)
    if [[ "$fit_ok" == False ]]; then
        echo "[WARN]  tumour prior of full/ was fitted below --prior-min-cov $TUMOR_PRIOR_MIN_COV (see tumor_atlas/full/tumor_prior.json)"
    fi

    # ---- 4: score against the tumour model (writes all five scores) ----
    for i in "${!samples[@]}"; do
        s=${samples[$i]}; o="$WORK/$s"
        model=full
        [[ "$BIASED" != 1 && "${categories[$i]}" == tumor_high ]] && model="without_$s"
        if [[ -f "$o/scores_tumor/read_scores.parquet" ]]; then
            echo "[SKIP]  $s tumour score"
        else
            rm -rf "$o/scores_tumor"
            echo "[START] $s tumour score (model $model, $TUMOR_MODEL)"
            tumor_args=(--tumor-atlas "$WORK/tumor_atlas/$model/tumor_atlas.feather")
            if [[ "$TUMOR_MODEL" == calibration ]]; then
                tumor_args+=(--tumor-calibration "$WORK/tumor_atlas/$model/tumor_calibration.json")
            else
                tumor_args+=(--tumor-prior auto --tumor-prior-min-cov "$TUMOR_PRIOR_MIN_COV")
            fi
            # shellcheck disable=SC2086
            nanoflux score -i "${calls[$i]}" --atlas "$ATLAS" -o "$o/scores_tumor" -c "${tumor_args[@]}" \
                --weight-score llr_tumor_control_per_call ${moments_arg[@]+"${moments_arg[@]}"} $SCORE_ARGS
        fi
    done
else
    echo "[WARN]  fewer than 2 tumor_high samples: no tumour model, the tumour scores are skipped"
fi

# ---- 5: evaluation ----
eval_sheet="$WORK/eval_sheet.tsv"
printf 'sample_id\tcategory\ttrue_label\treads_tsv\tscores\tscores_tumor\n' > "$eval_sheet"
for i in "${!samples[@]}"; do
    s=${samples[$i]}; o="$WORK/$s"
    tumor=""
    [[ -f "$o/scores_tumor/read_scores.parquet" ]] && tumor="$o/scores_tumor/read_scores.parquet"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "${categories[$i]}" "$(field $((i+1)) true_label)" \
        "$o/baseline/reads.tsv" "$o/scores/read_scores.parquet" "$tumor" >> "$eval_sheet"
done
echo "[START] evaluation"
"$PYTHON" "$HERE/evaluate_approaches.py" --sheet "$eval_sheet" --out "$WORK/results" --threads "$THREADS"
echo "[INFO] $(date)  finished. Table: $WORK/results/approaches_table.txt"
