process SC_RNASEQ_QC {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(counts)
    val(reference) // Path to sqlite of reference gff, generated with
    // ensembledb::ensDbFromGff
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: counts
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "QC", "h5ad")
    check = file("${meta.out}/${output}")
    extras = Utl.getName(module_number, meta, "QC_extras")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/sc_rnaseq_qc.log .
        cp -r ${meta.out}/${extras} .
        """
    } else {
        """
        mkdir ${extras}

        Rscript sc_rnaseq.R \\
          ${args} \\
          --thresholds_output ${extras}/mads_thresholds.tsv \\
          --plot TRUE \\
          --x_axis metric_x \\
          --y_axes sum,detected,subsets_mito_percent \\
          --plot_name ${extras}/diagnostics.png \\
          --loss_plot_name ${extras}/loss_diagnostics.png \\
          --discard \\
          --discard_out ${extras}/discarded_cells_output.tsv \\
          --input ${counts} \\
          --output ${output}

        get_nextflow_log.bash sc_rnaseq_qc.log
        """
    }
    //
}
