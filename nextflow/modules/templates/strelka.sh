configureStrelkaSomaticWorkflow.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    --indelCandidates !{manta_indels} \
    !{args} \
    --runDir !{out}

./runWorkflow.py

mv !{out}/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo $variant | sed 's/\.vcf\.gz//')

    vcf_info_add_tag -n SOURCE \\
        -d "Tool producing call" \\
        -b '.' \\
        -t String \\
        -a strelka2 \\
        -i $variant \\
        -o "!{module_number}-${base}_Strelka.vcf.gz"
done

cp .command.out strelka.log
