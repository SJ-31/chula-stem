#!/usr/bin/env bash
prefix="${1}"

for f in ../output/too_models/*/"${prefix}"-*hugo.csv; do
    dir=$(dirname "$f")
    dir="${dir}/bpformer"
    prediction="${dir}/bpformer-${prefix}.csv"
    report="${dir}/bpformer_report-${prefix}.txt"
    cm="${dir}/bpformer_confusion_matrix-${prefix}.csv"
    metrics="${dir}/bpformer_metrics-${prefix}.txt"
    ./predict_wrapper.bash -i "${f}" -p "${prediction}" -r "${report}" -c "${cm}" \
        -m "BPformer" -l "tumor_type" -e "${metrics}"
done
