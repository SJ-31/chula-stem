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
# Memory options/cores shouled be controlled by SLURM so don't need to be specified

mv !{out}/results/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    if [[ ! "$variant" == "candidateSmallIndels.vcf.gz" ]]; then
        base=$(echo $variant | sed 's/\.vcf\.gz//')
        name="!{module_number}-${base}_Strelka.vcf"

        rename_vcf.bash -v -i $variant -o tmp.vcf.gz \
            -n "${meta.RGSM_normal}" -t "${meta.RGSM_tumor}"

        vcf_info_add_tag -n SOURCE \
            -d "!{params.source_description}" \
            -b '.' \
            -t String \
            -a strelka2 \
            -i tmp.vcf.gz \
            -o "${name}"

        rm "$variant"
        bgzip "${name}"
    fi
done

get_nextflow_log.bash strelka.log
