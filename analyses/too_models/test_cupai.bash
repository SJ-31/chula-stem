#!/usr/bin/env bash
prefix="${1}"

for f in ../output/too_models/*/"$prefix"-*entrez.csv; do
    dir=$(dirname "$f")
    dir="${dir}/cupai"
    prediction="${dir}/cupai-${prefix}.csv"
    report="${dir}/cupai_report-${prefix}.txt"
    cm="${dir}/cupai_confusion_matrix-${prefix}.csv"
    metrics="${dir}/cupai_metrics-${prefix}.txt"

    ./predict_wrapper.bash -i "${f}" -p "${prediction}" -r "${report}" -c "${cm}" \
        -m "CUP-AI-Dx" -l "tumor_type" -e "${metrics}"
done
