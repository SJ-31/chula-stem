#!/usr/bin/env bash

while getopts "i:o:f:q:" f; do
    case "$f" in
        i) input=${OPTARG} ;;
        o) output=${OPTARG} ;;
        f) field=${OPTARG} ;;
        q) querystr=${OPTARG} ;;
        *)
            echo "Not supported"
            exit 1
            ;;
    esac
done

if [[ -z "$field" ]]; then
    field="ANN"
fi

header_str=$(bcftools head "${input}" | \
    grep "Consequence annotations from Ensembl VEP" | \
    sed -r 's/.*Format:(.*).*\">/\1/' | \
    sed 's/[|]/\t/g')

has_transcription_factors=$(echo "$header_str" | grep "TRANSCRIPTION_FACTORS")

if [[ -n "${has_transcription_factors}" ]]; then
    echo "Removing TRANSCRIPTION_FACTORS..."
    IFS=$'\t' read -r -a headers <<< "${header_str}"
    final="${#headers[@]}"
    exclude=$((final - 1))
    echo "${header_str}" | cut -f 1-"${exclude}" > "${output}"
    bcftools query "${input}" -f "%INFO/${field}" | \
        sed 's/[|]/\t/g' | \
        cut -f 1-"${exclude}" >> "${output}"
else
    echo "${header_str}" > "${output}"
    bcftools query "${input}" -f "%INFO/${field}" | \
        sed 's/[|]/\t/g' >> "${output}"
fi

if [[ -n "${querystr}" ]]; then
    mkdir temp
    echo -e "${querystr}" | sed 's/[][%]//g' | sed 's/[/]/_/g' > temp/query.tsv
    bcftools query "${input}" -f "${querystr}" >> temp/query.tsv
    mv "${output}" temp/to_join.tsv
    paste temp/query.tsv temp/to_join.tsv > "${output}"
    rm -R temp
fi
