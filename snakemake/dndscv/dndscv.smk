from pathlib import Path

cohort = config.get("cohort")
outdir = config.get("outdir", ".")

name = "dndscv" if not cohort else f"{cohort}-dndscv"
outdir = f"{outdir}/{name}"


rule dndscv:
    input:
        f"{outdir}/cohort_mutations.tsv",
    output:
        globaldnds=f"{outdir}/globaldnds.tsv",
        sel_cv=f"{outdir}/sel_cv.tsv",
        sel_loc=f"{outdir}/sel_loc.tsv",
        annotmuts=f"{outdir}/annotmuts.tsv",
        genemuts=f"{outdir}/genemuts.tsv",
        geneindels=f"{outdir}/geneindels.tsv",
        mle_submodel=f"{outdir}/mle_submodel.tsv",
        exclmuts=f"{outdir}/exclmuts.tsv",
        exclsamples=f"{outdir}/exclsamples.txt",
        nbreg=f"{outdir}/nbreg.rds",
        nbregind=f"{outdir}/nbreg.rds",
        poissmodel=f"{outdir}/poissmodel.rds",
        wrongmuts=f"{outdir}/wrongmuts.tsv",
    script:
        "dndscv_wrapper.R"


samples = {}
with open(config["sample_sheet"], "r") as f:
    lines = f.read().strip().splitlines()
    for line in lines:
        splits = line.split("\t")
        if len(splits) == 1:
            samples[Path(splits[0]).stem] = splits[0]
        else:
            samples[splits[0]] = splits[1]


rule get_mutations:
    output:
        f"{outdir}/cohort_mutations.tsv",
    run:
        for name, path in samples:
            query = f"{name}\t%CHROM\t%POS\t%REF\t%ALT"
            shell(f"bcftools query -f '{query}' {path} > {outdir}/{name}.tmp")
        shell(f"cat {outdir}/*tmp > {output[0]}")
