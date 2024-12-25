#!/usr/bin/env bash

rename=$(yq ".rename.grch38_nums" meta.yaml)
testdir=$(yq ".dir.pytest" meta.yaml)
a="${testdir}/vcf_test.vcf"
snpsift=$(yq ".tools.snpsift" meta.yaml)
java=$(yq ".java21" meta.yaml)
source="/data/project/stemcell/shannc/tests/wes_test/strelka2/4-somatic.snvs_Strelka.vcf.gz"
table="${testdir}/dndscv_mutants.tsv"

bcftools annotate --rename-chrs "$rename" "$source" > tmp.vcf
"$java" -jar "$snpsift" varType tmp.vcf > "$a"
bgzip -f "$a"
rm tmp.vcf
echo -e "chr\tpos\tref\talt" > "$table"
bcftools query \
    -i "SNP=1 || INS=1 || DEL=1" \
    -f "%CHROM\t%POS0\t%REF\t%ALT" \
    "$a.gz" >> "$table"

