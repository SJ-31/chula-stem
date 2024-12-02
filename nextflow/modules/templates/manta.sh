configManta.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    !{target_flag} \
    !{args} \
    --runDir !{out}

!{out}/runWorkflow.py

mv !{out}/results/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo $variant | sed -e 's/\.vcf\.gz//' -e 's/somatic\.//')
    name="!{prefix}-${base}_Manta.vcf"

    has_format=$( bcftools head "${variant}" | tail -n 1 | grep FORMAT )
    if [[ -n "${has_format}" && ! $variant =~ "diploid" ]]; then
        rename_vcf.bash -v -i $variant -n "!{meta.RGSM_normal}" -t "!{meta.RGSM_tumor}" | \
            vcf_info_add_tag.bash -n "!{params.source_name}" \
                -d "!{params.source_description}" \
                -b '.' \
                -t String \
                -a manta \
                -o "${name}"
    else
        vcf_info_add_tag.bash -n "!{params.source_name}" \
            -d "!{params.source_description}" \
            -b '.' \
            -t String \
            -a manta \
            -i "${variant}" \
            -o "${name}"
    fi
done

get_nextflow_log.bash manta.log
