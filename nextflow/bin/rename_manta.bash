#!/usr/bin/env bash

# If a letter is followed by a colon, the option expects an argument
while getopts "p:n:t:o:e:d:" f; do
    case "$f" in
        o) out=${OPTARG};;
        t) tumor=${OPTARG} ;;
        n) normal=${OPTARG} ;;
        e) source_name=${OPTARG} ;;
        d) source_description=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
        *)
            echo "Flag not recognized"
            exit 1
            ;;
    esac
done


mv "${out}"/results/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo "$variant" | sed -e 's/\.vcf\.gz//' -e 's/somatic\.//')
    name="${prefix}-${base}_Manta.vcf.gz"

    has_format=$( bcftools head "${variant}" | tail -n 1 | grep FORMAT )
    if [[ "${normal}" != "none" && -n "${has_format}" && ! $variant =~ "diploid" ]]; then
        rename_vcf.bash -v -i "$variant" -n "${normal}" -t "${tumor}" | \
            vcf_info_add_tag.bash -n "${source_name}" \
                -d "${source_description}" \
                -b '.' \
                -t String \
                -a manta \
                -o "${name}"
    else
        vcf_info_add_tag.bash -n "${source_name}" \
            -d "${source_description}" \
            -b '.' \
            -t String \
            -a manta \
            -i "${variant}" \
            -o "${name}"
    fi
done
