#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 5 || $# -gt 5 ]]; then
	echo "Usage: $0 <sample_id> <quantile> <min_cpg_filter> <short_read_filter> <metric>"
	echo "Example: $0 N2024-4049 0.6 50 2 nontumor_score_sum"
	exit 1
fi

sample_id=$1
quantile=$2
min_cpg_filter=$3
short_read_filter=$4
metric=$5

# Threshold min-CpG mode:
# - input: use the provided min_cpg_filter
# - zero: force threshold lookup from column 0
THRESHOLD_MIN_CPG_MODE="${THRESHOLD_MIN_CPG_MODE:-input}"
case "$THRESHOLD_MIN_CPG_MODE" in
	input|zero)
		;;
	*)
		echo "[ERROR] Unsupported THRESHOLD_MIN_CPG_MODE: $THRESHOLD_MIN_CPG_MODE (expected input or zero)" >&2
		exit 1
		;;
esac

lookup_quantile="$quantile"
lookup_min_cpg_filter="$min_cpg_filter"
if [[ "$THRESHOLD_MIN_CPG_MODE" == "zero" ]]; then
	lookup_min_cpg_filter="0"
fi

case "$metric" in
	nontumor_score_sum|nontumor_score_mean)
		;;
	*)
		echo "[ERROR] Unsupported metric: $metric (expected nontumor_score_sum or nontumor_score_mean)" >&2
		exit 1
		;;
esac

if [[ -z "${THRESHOLDS_CSV:-}" ]]; then
	THRESHOLDS_CSV="/home/lciernik/general_notebooks/thresholds_control_${metric}.csv"
fi

if [[ ! -f "$THRESHOLDS_CSV" ]]; then
	echo "[ERROR] Thresholds CSV not found: $THRESHOLDS_CSV" >&2
	exit 1
fi

score_threshold="$(python3 - "$THRESHOLDS_CSV" "$lookup_quantile" "$lookup_min_cpg_filter" "$metric" <<'PY'
import csv
import sys

csv_path, q_raw, cpg_raw, metric = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
q = float(q_raw)
cpg_col = str(int(cpg_raw))

with open(csv_path, newline="") as f:
	reader = csv.DictReader(f)
	if cpg_col not in reader.fieldnames and "0" in reader.fieldnames:
		cpg_col = "0"
	if cpg_col not in reader.fieldnames:
		raise SystemExit(f"Missing column '{cpg_col}' in {csv_path}")

	best = None
	best_diff = None
	for row in reader:
		try:
			q_row = float(row["quantile"])
		except (KeyError, ValueError):
			continue

		diff = abs(q_row - q)
		if best is None or diff < best_diff:
			best = row
			best_diff = diff

	if best is None or best_diff is None or best_diff > 1e-9:
		raise SystemExit(f"No quantile match for {q_raw} in {csv_path}")

	val = best.get(cpg_col, "")
	if val == "":
		raise SystemExit(
			f"Empty threshold for quantile={q_raw}, min_cpg_filter={cpg_col} in {csv_path}"
		)

	print(val)
PY
)"

setting_id="metric_${metric}_quantile_${quantile}_min_cpgs_${min_cpg_filter}_short_read_filter_min_covered_cpgs_${short_read_filter}_v2"
bam_path="/home/lciernik/shared_space_nanopore/calls/${sample_id}/aligned_to_13v2.bam"
read_filter_feather="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/aligned_to_hg38.nontumor_pooled_ref_scored.feather"
sample_out_dir="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/${setting_id}"
if [[ "$THRESHOLD_MIN_CPG_MODE" == "zero" ]]; then
	sample_out_dir="${sample_out_dir}_quantile_from_total_distribution"
fi
mkdir -p "$sample_out_dir"

echo "[INFO] Using score threshold ${score_threshold} from metric=${metric}, input_quantile=${quantile}, lookup_quantile=${lookup_quantile}, input_min_cpg_filter=${min_cpg_filter}, lookup_min_cpg_filter=${lookup_min_cpg_filter}"

if [[ ! -f "$bam_path" ]]; then
    echo "[ERROR] BAM file not found: $bam_path" >&2
    exit 1
fi

if [[ ! -f "$read_filter_feather" ]]; then
    echo "[ERROR] Read-filter feather file not found: $read_filter_feather" >&2
    exit 1
fi

# If predictions are already present, skip processing
if [[ -f "${sample_out_dir}/predicted_labels.json" ]]; then
	echo "[INFO] Predictions already exist for sample $sample_id with setting $setting_id, skipping processing"
	exit 0
fi

tmp_sample_dir="$(mktemp -d "$sample_out_dir/.nanoflux_prepare.${sample_id}.XXXXXX")"
tmp_bam="$tmp_sample_dir/aligned_to_13v2.bam"

echo "[INFO] Copying BAM to temporary output workspace: $tmp_bam"
if ! cp --reflink=auto "$bam_path" "$tmp_bam" 2>/dev/null; then
    cp "$bam_path" "$tmp_bam"
fi

(
    trap 'rm -rf "$tmp_sample_dir"' EXIT

	nanoflux prepare \
		-i "$tmp_bam" \
		-o "$sample_out_dir" \
		-c \
		--skip-alignment \
		--read-filter-feather "$read_filter_feather" \
		--non-tumor-score-threshold "$score_threshold" \
		--score-filter-min-cpgs "$min_cpg_filter" \
		--short-read-min-covered-cpgs "$short_read_filter" \
		--feather-cols read_name "$metric" covered_cpgs \
		-t 32
)
nanoflux infer -i "${sample_out_dir}/methylation.bed" -o "$sample_out_dir"

# remove the methylation.bed file to save space, it can be regenerated from the BAM if needed
rm -f "${sample_out_dir}/methylation.bed"