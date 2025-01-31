#!/usr/bin/env bash

while getopts "l:e:i:p:c:r:m:" f; do
    case "$f" in
        i) input=${OPTARG} ;; # Stores the argument
        p) prediction=${OPTARG} ;;
        r) report=${OPTARG} ;;
        m) model=${OPTARG} ;;
        c) cm=${OPTARG} ;;
        e) metrics=${OPTARG} ;;
        l) label=${OPTARG} ;;
        *) echo "Invalid flag"; exit 1;;
    esac
done

input=$(realpath "${input}")
prediction=$(realpath "${prediction}")
report=$(realpath "${report}")
cm=$(realpath "${cm}")
metrics=$(realpath "${metrics}")

case "${model}" in
    "CUP-AI-Dx") dir="../../extras/CUP-AI-Dx";;
    "BPformer") dir="../../extras/BPformer";;
    *) echo "Invalid model"; exit 1;;
esac

pushd "${dir}" || echo "No directory!"
python run.py -i "${input}" -o "${prediction}" -r "${report}"  -c "${cm}" -l "${label}" \
    -m "${metrics}"
popd
