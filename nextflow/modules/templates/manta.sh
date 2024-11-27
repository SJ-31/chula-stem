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
    name="!{module_number}-!{meta.filename}-${base}_Manta.vcf"

    vcf_info_add_tag -n "!{params.source_name}" \
        -d "!{params.source_description}" \
        -b '.' \
        -t String \
        -a manta \
        -i $variant \
        -o "${name}"

    rm "$variant"
    bgzip "${name}"
done

get_nextflow_log.bash manta.log
