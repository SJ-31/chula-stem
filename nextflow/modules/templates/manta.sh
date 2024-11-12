bcftools index !{normal}
bcftools index !{tumor}

configManta.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
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
        -a manta \\
        -i $variant \\
        -o "!{module_number}-${base}_Manta.vcf.gz"
done

bcftools annotate temp.vcf.gz TODO
cp .command.out manta.log
