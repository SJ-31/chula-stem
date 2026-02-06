import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell
def _():
    from yte import process_yaml

    import anndata as ad
    import marimo as mo
    from pathlib import Path
    import scanpy as sc
    import plotnine as gg
    from chula_stem.sc_rnaseq import pca_to_leiden

    with open("./cellranger_config.yaml", "r") as f:
        env = process_yaml(f)
        print(env)

    cohort_dir = Path(env["outdir"]) / "cohort"
    plot_out = cohort_dir / "plots"
    return ad, env, mo, plot_out, sc


@app.cell
def _():
    import functions as fn
    return (fn,)


@app.cell
def _(ad, env, fn):
    combined: ad.AnnData = fn.prepare_data(env)
    return (combined,)


@app.cell
def _(combined: "ad.AnnData"):
    print(combined.obs.columns)
    combined.obsm["mads"]
    return


@app.cell
def _(combined: "ad.AnnData", sc):
    sc.pp.subsample(combined, copy=True, fraction=0.1).write_h5ad("./test.h5ad")
    combined.var
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Quality control
    """)
    return


@app.cell
def _(combined: "ad.AnnData", fn, plot_out):
    (plot_out / "qc").mkdir(exist_ok=True, parents=True)
    for patient in combined.obs["patient"].unique():
        if (plot_out / "qc" / f"{patient}.pdf").exists():
            continue
        qc_mads, mads_compare = fn.qc_plot_patient(combined.to_memory(), patient)
        (qc_mads / mads_compare).save(plot_out / "qc" / f"{patient}.pdf")
    return


@app.cell
def _(combined: "ad.AnnData", env, fn):
    filtered, failed = fn.mads_filter_outliers(
        combined.to_memory(), filters=env["mads_thresholds"], reduction="all"
    )
    print("Shape of filtered: {}".format(filtered))
    print("Shape of failed: {}".format(failed))

    filtered.write_h5ad(env["files"]["passed_qc"])
    failed.write_h5ad(env["files"]["failed_qc"])
    return


if __name__ == "__main__":
    app.run()
