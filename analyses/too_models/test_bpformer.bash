#!/usr/bin/env bash

for f in ../output/too_models/*/*hugo.csv; do
    dir=$(dirname "$f")
    prediction="${dir}/bpformer.csv"
    report="${dir}/bpformer_report.txt"
    cm="${dir}/bpformer_confusion_matrix.csv"
    ./predict_wrapper.bash -i "${f}" -p "${prediction}" -r "${report}" -c "${cm}" \
        -m "BPformer" -l "tumor_type"
done
