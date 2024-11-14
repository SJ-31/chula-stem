configManta.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    !{args} \
    --runDir !{out}

!{out}/runWorkflow.py

mv !{out}/results/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo $variant | sed 's/\.vcf\.gz//')
    name="!{module_number}-${base}_Manta.vcf"

    vcf_info_add_tag -n SOURCE \
        -d "!{params.source_description}" \
        -b '.' \
        -t String \
        -a manta \
        -i $variant \
        -o "${name}"

    bgzip "${name}"
done

cp .command.out manta.log
