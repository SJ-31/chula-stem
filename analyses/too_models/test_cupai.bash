#!/usr/bin/env bash

for f in ../output/too_models/*/*entrez.csv; do
    dir=$(dirname "$f")
    prediction="${dir}/cupai.csv"
    report="${dir}/cupai_report.txt"
    cm="${dir}/cupai_confusion_matrix.csv"
    metrics="${dir}/cupai_metrics.txt"
    ./predict_wrapper.bash -i "${f}" -p "${prediction}" -r "${report}" -c "${cm}" \
        -m "CUP-AI-Dx" -l "tumor_type" -e "${metrics}"
done
