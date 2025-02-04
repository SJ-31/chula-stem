#!/usr/bin/env bash

while getopts "o:l:e:i:p:c:r:m:" f; do
    case "$f" in
        i) input=$(realpath "${OPTARG}") ;;
        p) prediction=$(realpath "${OPTARG}") ;;
        r) report=$(realpath "${OPTARG}") ;;
        m) model=${OPTARG} ;;
        c) cm=$(realpath "${OPTARG}") ;;
        e) metrics=$(realpath "${OPTARG}") ;;
        l) label=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
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
    python run.py -i "${input}" -o "${prediction}" -r "${report}" \
        -c "${cm}" -l "${label}" \
        -m "${metrics}"
fi
popd
