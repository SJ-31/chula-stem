// All-in-one module for QC on SC-RNASeq data
// Carries out both doublet detection and filtering for low-quality cells,
//      producing diagnostic plots

process SC_RNASEQ_QC {
    ext version: "1"

    input:
    tuple val(meta), path(counts)
    val reference // Path to sqlite of reference gff, generated with

    output:
    // ensembledb::ensDbFromGff

    //
    tuple val(meta), path(output), emit: counts
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    output = Utl.getName(meta, "QC", "h5ad")
    extras = Utl.getName(meta, "QC_extras")
    args = task.ext.args.join(" ")

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
