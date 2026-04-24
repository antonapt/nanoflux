#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
	echo "Usage: $0 <sample_id> <quantile> <min_cpg_filter> <metric>"
	echo "Example: $0 N2024-4049 0.6 50 nontumor_score_sum"
	exit 1
fi

sample_id=$1
quantile=$2
min_cpg_filter=$3
metric=$4

case "$metric" in
	nontumor_score_sum|nontumor_score_mean) ;;
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

setting_id="metric_${metric}_quantile_${quantile}_min_cpgs_${min_cpg_filter}_v2"
bam_path="/home/lciernik/shared_space_nanopore/calls/${sample_id}/aligned_to_13v2.bam"
read_filter_feather="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/aligned_to_hg38.nontumor_pooled_ref_scored.feather"
sample_out_dir="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/${setting_id}"
mkdir -p "$sample_out_dir"

echo "[INFO] Running approach 2 for sample=$sample_id metric=$metric quantile=$quantile min_cpg_filter=$min_cpg_filter"

if [[ ! -f "$bam_path" ]]; then
	echo "[ERROR] BAM file not found: $bam_path" >&2
	exit 1
fi

if [[ ! -f "$read_filter_feather" ]]; then
	echo "[ERROR] Read-filter feather file not found: $read_filter_feather" >&2
	exit 1
fi

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
		--thresholds-file "$THRESHOLDS_CSV" \
		--quantile "$quantile" \
		--score-filter-min-cpgs "$min_cpg_filter" \
		--feather-cols read_name "$metric" covered_cpgs \
		-t 32
)
nanoflux infer -i "${sample_out_dir}/methylation.bed" -o "$sample_out_dir"

rm -f "${sample_out_dir}/methylation.bed"
