#!/usr/bin/env bash
prefix="${1}"

for f in ../output/too_models/*/$prefix*entrez.csv; do
    dir=$(dirname "$f")
    dir="${dir}/cupai"
    for m in "inception" "cnn" "resnet"; do
        p="cupai-${prefix}-${m}"
        ./predict_wrapper.bash -i "${f}" -m "CUP-AI-Dx" \
            -l "tumor_type" -s "${m}" -p "${p}-" -o "${dir}"
    done
done
