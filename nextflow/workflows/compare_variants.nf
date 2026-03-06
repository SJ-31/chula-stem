include { VCFEVAL } from "../modules/vcfeval.nf"

workflow compare_variants {

    main:
    input = channel.fromPath(params.manifest)
        .splitCsv(header: true)
        .map { it -> [["id": it.sample,
                "out": "${params.outdir}/${it.sample}",
                "type": it.source,
                "log": "${params.logdir}/${it.sample}/${it.source}",
                "filename": "${it.sample}_${it.source}",
            ], [file(it.old_calls), file(it.new_calls)]] }

    VCFEVAL(input, params.ref.genome_sdf, params.ref.targets, 1)
}
