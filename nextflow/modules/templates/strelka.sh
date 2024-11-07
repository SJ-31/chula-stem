configureStrelkaSomaticWorkflow.py \
    --normalBam !{normal} \
    --tumorBam !{tumor} \
    --referenceFasta !{reference} \
    --indelCandidates !{manta_indels} \
    !{exome_flag} \
    --runDir !{out}

./runWorkflow.py

mv !{out}/variants/*.vcf.gz .
for variant in *.vcf.gz; do
    base=$(echo $variant | sed 's/\.vcf\.gz//')
    mv $variant "!{module_number}-${base}_Strelka.vcf.gz"
done

bcftools annotate temp.vcf.gz TODO
cp .command.out strelka.log
