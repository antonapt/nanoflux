#!/usr/bin/env bash
# Score every sample against the CSF cfDNA atlas. This is exactly what was run
# on the nine samples s_6 ... s_21 (see the comments for what happens inside).
#
# Usage:
#   scripts/score_samples.sh [samples_dir] [atlas.feather] [out_dir]
#
# Defaults are the paths used on the analysis Mac:
#   samples_dir  /Users/lciernik/Documents/TUB/projects/dna_methylation/read_filtering/data/samples
#   atlas        /Users/lciernik/Documents/TUB/projects/dna_methylation/read_filtering/data/reference_csf/reference.feather
#   out_dir      <samples_dir>            (scores land in <samples_dir>/<sample>/scores/)
#
# Input per sample:  <samples_dir>/<sample>/extracted.tsv
#   The modkit read-calls table of the sample aligned to hg38 (one row per read
#   and CpG: read_id, ref_position, chrom, ref_strand, call_prob, call_code, ...).
#   It must be on hg38 because the atlas is on hg38. Read names are the join key
#   to the chm13v2 pipeline later, so the alignment used for scoring does not
#   have to be the one used for prepare.
#
# What `nanoflux score` does with it, step by step:
#   1. Loads the atlas: one row per CpG with chrom_x, start_genomic (1-based
#      plus-strand C), methylated_pooled and total_pooled.
#   2. Fits the Beta prior to the well-covered atlas sites (>= 100 pooled reads):
#      p = (methylated + alpha) / (total + alpha + beta). On this atlas
#      alpha = 1.45, beta = 0.68, so shallow sites are pulled toward 68 percent
#      methylated instead of toward 50 percent. (Override with --prior A B.)
#   3. Reads extracted.tsv in chunks, keeps primary alignments and the calls
#      "m" (5mC) and "-" (unmethylated); "h" (5hmC) is dropped. Each call is
#      matched to the atlas: plus-strand reads at ref_position + 1, minus-strand
#      reads at ref_position (modkit reports the G there). About 99.7 percent
#      of calls match.
#   4. Fits the caller moments on the sample itself: the mean and variance of
#      the read's probability of methylation at truly methylated and truly
#      unmethylated sites (mu1, v1, mu0, v0), by two linear regressions of
#      r and r^2 on p. Written to moments.json. (Reuse across samples with
#      --moments <file>.)
#   5. Scores every read:
#        obs = sum over its CpGs of  log(1-p) + r * logit(p)
#        exp = the same with the expected r for a read that follows the atlas
#        var = the variance of that sum
#        z   = (obs - exp) / sqrt(var)
#      z near 0: the read looks like the atlas; strongly negative: it does not.
#   6. Turns z into a weight without any cutoff:
#        weight = sigmoid(-z / temperature)         (temperature 1 by default)
#      so a read at z = 0 gets 0.5, at z = -2 about 0.88, at z = -4 about 0.98.
#      Only the ratio of weights between reads at the same CpG matters later.
#
# Output per sample (in <out_dir>/<sample>/scores/):
#   read_scores.parquet   read_id, n_calls, n_cpg, n_bin, z, llr_per_call, weight
#   moments.json          the prior and the fitted caller moments
#   score_summary.json    counts and per-bin z quantiles and mean weights

set -euo pipefail

SAMPLES_DIR=${1:-/Users/lciernik/Documents/TUB/projects/dna_methylation/read_filtering/data/samples}
ATLAS=${2:-/Users/lciernik/Documents/TUB/projects/dna_methylation/read_filtering/data/reference_csf/reference.feather}
OUT_DIR=${3:-$SAMPLES_DIR}

command -v nanoflux >/dev/null || { echo "nanoflux not on PATH; activate the environment" >&2; exit 1; }

for dir in "$SAMPLES_DIR"/*/; do
    sample=$(basename "$dir")
    calls="$dir/extracted.tsv"
    out="$OUT_DIR/$sample/scores"
    if [[ ! -f "$calls" ]]; then
        echo "[SKIP] $sample: no extracted.tsv"
        continue
    fi
    if [[ -f "$out/read_scores.parquet" ]]; then
        echo "[SKIP] $sample: $out/read_scores.parquet exists"
        continue
    fi
    echo "[RUN]  $sample"
    nanoflux score -i "$calls" --atlas "$ATLAS" -o "$out" -c
done
