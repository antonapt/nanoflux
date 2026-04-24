#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 || $# -gt 2 ]]; then
	echo "Usage: $0 <sample_id> <short_read_min_covered_cpgs>"
	echo "Example: $0 N2024-4049 50"
	exit 1
fi

sample_id=$1
short_read_min_covered_cpgs=$2

bam_path="/home/lciernik/shared_space_nanopore/calls/${sample_id}/aligned_to_13v2.bam"
read_filter_feather="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/aligned_to_hg38.nontumor_pooled_ref_scored.feather"
setting_id="short_read_filter_min_cpgs_${short_read_min_covered_cpgs}"
sample_out_dir="/home/lciernik/storage/data/wgbs_pp_data/${sample_id}/${setting_id}"
mkdir -p "$sample_out_dir"

echo "[INFO] Filtering reads with minimum covered CpGs: ${short_read_min_covered_cpgs}"

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
		--short-read-min-covered-cpgs "$short_read_min_covered_cpgs" \
		-t 32
)
nanoflux infer -i "${sample_out_dir}/methylation.bed" -o "$sample_out_dir"

# remove the methylation.bed file to save space, it can be regenerated from the BAM if needed
rm -f "${sample_out_dir}/methylation.bed"