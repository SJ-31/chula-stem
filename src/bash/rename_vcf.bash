#!/usr/bin/bash
#
# Rename samples in the provided paired-normal vcf file to NORMAL, TUMOR and reorder
# the vcf to be NORMAL\tTUMOR
#
reverse=0
while getopts "vi:o:n:t:" f; do
    case "$f" in
        i)
            input="${OPTARG}"
            ;;
        o)
            output="${OPTARG}"
            ;;
        n)
            normal_name="${OPTARG}"
            ;;
        t)
            tumor_name="${OPTARG}"
            ;;
        v)
            reverse=1
            ;;
        *)
            echo "Unknown flag"
            exit 1
    esac
done

if [[ ${reverse} == 0 ]]; then
    rename_str="${normal_name} NORMAL\n${tumor_name} TUMOR"
    order_str="NORMAL,TUMOR"
else
    rename_str="NORMAL ${normal_name}\nTUMOR ${tumor_name}"
    order_str="${normal_name},${tumor_name}"
fi

echo -e "${rename_str}" > rename.txt
if [[ -z "${output}" ]]; then
    bcftools reheader -s rename.txt "${input}" | bcftools view -s "${order_str}"
else
    bcftools reheader -s rename.txt "${input}" | \
        bcftools view -s "${order_str}" -O z > "${output}"
fi
