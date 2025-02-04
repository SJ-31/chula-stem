#!/usr/bin/env bash
prefix="${1}"


for f in ../output/too_models/*/"$prefix"-*tulip.csv; do
    dir=$(dirname "$f")
    dir="${dir}/tulip"
    ./predict_wrapper.bash -i "${f}" -o "${dir}" -m "TULIP"
done
