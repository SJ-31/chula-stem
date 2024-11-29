#!/usr/bin/env bash

# If a letter is followed by a colon, the option expects an argument
#

while getopts "i:o:b:n:d:t:a:" f; do
    case "$f" in
        i) input=${OPTARG} ;;
        o) output=${OPTARG} ;;
        b) number=${OPTARG} ;;
        n) name=${OPTARG} ;;
        d) description=${OPTARG} ;;
        t) type=${OPTARG} ;;
        a) default=${OPTARG} ;;
        *) echo "Flag not provided"
           exit 1 ;;
    esac
done

echo "##INFO=<ID=${name},Number=${number},Type=${type},Description=\"${description}\">" > line.txt

if [[ -z "${input}" ]]; then
    input=$(cat - )
else
    input=$(bcftools view "${input}")
fi

echo "$input" | bcftools annotate -h line.txt | \
    awk -v tn="${name}" -v val="${default}" \
    'BEGIN {OFS="\t"} {if ($1 !~ "^#") {$8=$8";"tn"="val} {print}}' | \
    bcftools view -O z > "${output}"

## Output is automatically bgzipped

rm line.txt
