#!/usr/bin/env bash

while getopts "i:p:c:r:m:" f; do
    case "$f" in
        i) input=${OPTARG} ;; # Stores the argument
        p) prediction=${OPTARG} ;;
        r) report=${OPTARG} ;;
        m) model=${OPTARG} ;;
        c) cm=${OPTARG} ;;
        l) label=${OPTARG} ;;
        *) echo "Invalid flag"; exit 1;;
    esac
done

input=$(realpath "${input}")
prediction=$(realpath "${prediction}")
report=$(realpath "${report}")
cm=$(realpath "${cm}")

case "${model}" in
    "CUP-AI-Dx") dir="../../extras/CUP-AI-Dx";;
    "BPformer") dir="../../extras/BPformer";;
    *) echo "Invalid model"; exit 1;;
esac

pushd "${dir}" || echo "No directory!"
python run.py -i "${input}" -o "${prediction}" -r "${report}"  -c "${cm}" -l "${label}"
popd
