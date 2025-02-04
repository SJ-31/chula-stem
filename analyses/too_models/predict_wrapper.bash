#!/usr/bin/env bash

while getopts "o:l:i:s:p:m:" f; do
    case "$f" in
        i) input=$(readlink -m "${OPTARG}") ;;
        m) model=${OPTARG} ;;
        s) model_spec=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
        l) label=${OPTARG} ;;
        o) outdir=$(readlink -m "${OPTARG}") ;;
        *) echo "Invalid flag"; exit 1;;
    esac
done

case "${model}" in
    "CUP-AI-Dx") dir="../../extras/CUP-AI-Dx";;
    "BPformer") dir="../../extras/BPformer";;
    "TULIP") dir="../../extras/TULIP";;
    *) echo "Invalid model"; exit 1;;
esac

pushd "${dir}" || echo "No directory!"
if [[ "${model}" == "TULIP" ]]; then
    python tulip.py -i "${input}" -t 17 -o "${outdir}"
else
    python run.py -i "${input}" -o "${outdir}" -p "${prefix}" \
        -l "${label}" --model "${model_spec}"
fi
popd
