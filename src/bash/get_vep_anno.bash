#!/usr/bin/env bash

input="$1"
output="$2"
if [[ -z "$3" ]]; then
    field="ANN"
else
    field="$3"
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
