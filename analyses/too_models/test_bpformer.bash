#!/usr/bin/env bash
prefix="${1}"

for f in ../output/too_models/*/"${prefix}"*hugo.csv; do
    dir=$(dirname "$f")
    dir="${dir}/bpformer"
    ./predict_wrapper.bash -i "${f}" \
        -m "BPformer" -l "tumor_type" -o "${dir}" -p "bpformer-${prefix}-"
done
