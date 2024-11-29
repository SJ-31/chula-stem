#!/usr/bin/env bash

cnvkit.py batch \
        !{tumor} \
    -r !{cnn_reference} \
    -d !{out}

if [[ !{clonal} == "true" ]]; then
    cnvkit.py call \
        !{out}/!{out}.call.cns \
        --vcf !{snps} \
        --purity !{purity_val} \
        --ploidy !{ploidy_val} \
        --method clonal \
        --output !{out}.call.clonal.cns
    mv !{out}.call.clonal.cns !{out}
fi

cd !{out}
for file in !{tumor.baseName}*; do
    prefix=$(echo "${file}" | sed 's/!{tumor.baseName}//')
    mv "${file}" "!{out}${prefix}"
done
cd ..

get_nextflow_log.bash cnvkit.log
