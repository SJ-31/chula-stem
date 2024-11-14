tabix -f !{manta_indels}

configureStrelkaSomaticWorkflow.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    --indelCandidates !{manta_indels} \
    !{args} \
    --runDir !{out}

!{out}/runWorkflow.py -m local
# Memory options/cores shouled be controlled by SLURM so don't need to be specified

mv !{out}/results/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo $variant | sed 's/\.vcf\.gz//')
    name="!{module_number}-${base}_Strelka.vcf"

    vcf_info_add_tag -n SOURCE \
        -d "!{params.source_description}" \
        -b '.' \
        -t String \
        -a strelka2 \
        -i $variant \
        -o "${name}"

    bgzip "${name}"
done

cp .command.out strelka.log
