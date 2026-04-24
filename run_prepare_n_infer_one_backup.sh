#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 || $# -gt 3 ]]; then
	echo "Usage: $0 <sample_id> <score_threshold> <min_cpg_filter>"
	echo "Example: $0 N2024-4049 -0.78 50"
	exit 1
fi

sample_id=$1
score_threshold=$2
min_cpg_filter=$3

setting_id="score_${score_threshold}_min_cpgs_${min_cpg_filter}"
bam_path="/home/lciernik/shared_space_nanopore/calls/${sample_id}/aligned_to_13v2.bam"
read_filter_feather="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/aligned_to_hg38.nontumor_pooled_ref_scored.feather"
sample_out_dir="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/${setting_id}"
mkdir -p "$sample_out_dir"

if [[ ! -f "$bam_path" ]]; then
    echo "[ERROR] BAM file not found: $bam_path" >&2
    exit 1
fi

if [[ ! -f "$read_filter_feather" ]]; then
    echo "[ERROR] Read-filter feather file not found: $read_filter_feather" >&2
    exit 1
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
		--min-covered-cpgs "$min_cpg_filter" \
		-t 32
)
nanoflux infer -i "${sample_out_dir}/methylation.bed" -o "$sample_out_dir"