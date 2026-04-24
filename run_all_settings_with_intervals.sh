#!/usr/bin/env bash

set -euo pipefail

all_samples=(
    N2023-3932
    N2023-3951
    N2023-4036
    N2023-4287
    N2025-3860
    N2025-3890
    N2025-4007
    N2025-4009
    N2025-4013
    N2025-4032
    N2025-4033
    N2025-4082
    N2025-4126
    N2021-1045
    N2023-0644
    N2025-0190
    N2025-0749
    N2025-3258
    N2025-4036
    N2020-0267
    N2025-2532
    N2019-2441
    N2024-4045
    N2024-4049
    N2020-2407
    N2023-0800
    N2024-1640
    N2023-3676
    N2025-3202
    N2024-2748
    N2022-0550
    N2024-4431
    N2021-0414
    N2023-1611
    N2025-2614
    A2023-0050
    N2025-4241
    N2024-0965
    N2025-3850
)
all_quantiles=( 0.05 0.1 0.2 0.4 0.7 0.8 0.975 )
all_metrics=('nontumor_score_mean' 'nontumor_score_sum')


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="$SCRIPT_DIR/run_prepare_n_infer_with_intervals_one_setting.sh"

if [[ ! -x "$RUNNER" ]]; then
    echo "Error: runner script is missing or not executable: $RUNNER" >&2
    exit 1
fi

# Tune parallelism with MAX_JOBS env var, e.g. MAX_JOBS=6 ./run_all_settings_with_intervals.sh
MAX_JOBS=6
MAX_JOBS="${MAX_JOBS:-$(nproc)}"
if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || [[ "$MAX_JOBS" -lt 1 ]]; then
    echo "Error: MAX_JOBS must be a positive integer, got '$MAX_JOBS'" >&2
    exit 1
fi

echo "[INFO] Running all combinations with MAX_JOBS=$MAX_JOBS"

run_stamp="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$SCRIPT_DIR/logs/run_${run_stamp}"
STATUS_DIR="$LOG_DIR/status"
mkdir -p "$LOG_DIR" "$STATUS_DIR"
export LOG_DIR STATUS_DIR
echo "[INFO] Writing per-job logs to: $LOG_DIR"

job_count=$(( ${#all_samples[@]} * ${#all_quantiles[@]} * ${#all_metrics[@]} ))

for sample_id in "${all_samples[@]}"; do
    for quantile in "${all_quantiles[@]}"; do
        for metric in "${all_metrics[@]}"; do
            printf '%s\t%s\t%s\n' "$sample_id" "$quantile" "$metric"
        done
    done
done |
xargs -P "$MAX_JOBS" -n 3 bash -c '
    set -uo pipefail
    sample_id="$1"
    quantile="$2"
    metric="$3"
    job_key="${sample_id}__q${quantile}__${metric}"
    log_file="$LOG_DIR/${sample_id}/${job_key}.log"
    status_file="$STATUS_DIR/${job_key}.status"
    mkdir -p "$LOG_DIR/${sample_id}"

    set +e
    echo "[START] sample=$sample_id quantile=$quantile metric=$metric" >"$log_file"
    "$0" "$sample_id" "$quantile" "$metric" >>"$log_file" 2>&1
    exit_code=$?

    if [[ "$exit_code" -eq 0 ]]; then
        echo "[DONE]  sample=$sample_id quantile=$quantile metric=$metric" >> "$log_file"
        echo "ok" > "$status_file"
    else
        echo "[FAILED] sample=$sample_id quantile=$quantile metric=$metric exit_code=$exit_code" >> "$log_file"
        echo "failed" > "$status_file"
    fi

    exit 0
' "$RUNNER"

failed_jobs=$(
    find "$STATUS_DIR" -type f -name '*.status' -print0 2>/dev/null \
        | xargs -0 -r cat 2>/dev/null \
        | awk '$0=="failed" {c++} END {print c+0}'
)
successful_jobs=$((job_count - failed_jobs))

echo "[INFO] Finished $job_count jobs"
echo "[INFO] Successful jobs: $successful_jobs"
echo "[INFO] Failed jobs: $failed_jobs"

if [[ "$failed_jobs" -gt 0 ]]; then
    echo "[INFO] Failed job logs:"
    grep -R -l '^failed$' "$STATUS_DIR" 2>/dev/null \
        | sed 's#.*/##; s#\.status$#.log#' \
        | while IFS= read -r f; do echo "  $LOG_DIR/$f"; done
    exit 1
fi
