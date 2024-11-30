tabix -f !{manta_indels}

configureStrelkaSomaticWorkflow.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    --indelCandidates !{manta_indels} \
    !{target_flag} \
    !{args} \
    --runDir !{out}

!{out}/runWorkflow.py -m local
# Memory options/cores should be controlled by SLURM so don't need to be specified

mv !{out}/results/variants/*.vcf.gz .
for variant in somatic*.vcf.gz; do
    base=$(echo $variant | sed -e 's/\.vcf\.gz//' -e 's/somatic\.//')
    name="!{prefix}-${base}_Strelka.vcf"

    rename_vcf.bash -v -i $variant -o tmp.vcf.gz \
        -n "!{meta.RGSM_normal}" -t "!{meta.RGSM_tumor}"

    vcf_info_add_tag.bash -n "!{params.source_name}" \
        -d "!{params.source_description}" \
        -b '.' \
        -t String \
        -a strelka2 \
        -i tmp.vcf.gz \
        -o "${name}"

    rm "$variant"
    bgzip "${name}"
done

get_nextflow_log.bash strelka.log
