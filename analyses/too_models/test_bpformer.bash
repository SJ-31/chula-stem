#!/usr/bin/env bash

for f in ../output/too_models/*/*hugo.csv; do
    dir=$(dirname "$f")
    prediction="${dir}/bpformer.csv"
    report="${dir}/bpformer_report.txt"
    cm="${dir}/bpformer_confusion_matrix.csv"
    metrics="${dir}/bpformer_metrics.txt"
    ./predict_wrapper.bash -i "${f}" -p "${prediction}" -r "${report}" -c "${cm}" \
        -m "BPformer" -l "tumor_type" -e "${metrics}"
done
