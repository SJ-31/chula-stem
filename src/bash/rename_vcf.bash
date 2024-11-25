#!/usr/bin/bash
#
# Rename samples in the provided paired-normal vcf file to NORMAL, TUMOR and reorder
# the vcf to be NORMAL\tTUMOR
#
while getopts "i:o:n:t:" f; do
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
        *)
            echo "Unknown flag"
            exit 1
    esac
done

echo -e "${normal_name} NORMAL\n${tumor_name} TUMOR" > rename.txt
bcftools reheader -s rename.txt "${input}" | \
    bcftools view -s NORMAL,TUMOR -O z > "${output}"
