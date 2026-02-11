import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell
def _():
    import sys
    from pathlib import Path

    import anndata as ad
    import marimo as mo
    import plotnine as gg
    import scanpy as sc
    from yte import process_yaml

    return Path, ad, mo, process_yaml, sys


@app.cell
def _(Path, mo):
    if Path.home() == "/home/shannc":
        is_test = mo.ui.dropdown(
            options=[True, False], value=True, label="Use test data?"
        )
        is_test
    else:
        is_test = None
    return (is_test,)


@app.cell(hide_code=True)
def _(Path, is_test, process_yaml, sys):
    use_test_data = is_test.value if is_test is not None else False
    if use_test_data:
        sys.argv.append("test=True")
    with open("./cellranger_config.yaml", "r") as f:
        env = process_yaml(f)
        print(env)

    plot_out = Path(env["outdir"]) / "plots"
    return env, plot_out


@app.cell
def _():
    import functions as fn

    return (fn,)


@app.cell
def _(ad, env, fn):
    combined: ad.AnnData = fn.data_import(env)
    return (combined,)


@app.cell
def _(combined: "ad.AnnData"):
    print(combined.obs.columns)
    combined.obsm["mads"]
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Quality control
    """)
    return


@app.cell
def _(ad, env, fn):
    adata = ad.read_h5ad(env["files"]["passed_qc"])
    fn.add_saved_dr(adata, env)
    return (adata,)


@app.cell
def _(mo):
    mo.md(r"""
    # Dimensionality reduction (unintegrated)
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## T-SNE
    """)
    return


@app.cell
def _(adata, fn, plot_out):
    # TODO: add this to the main notebook
    tsne_slider, display_tsne = fn.make_dr_slider(
        adata,
        "t-sne",
        plot_out / "tsne_unintegrated",
        ["patient", "type"],
        ext="svg",
        theme_kws={"figure_size": (10, 5)},
    )
    tsne_slider
    return display_tsne, tsne_slider


@app.cell
def _(display_tsne, tsne_slider):
    display_tsne(tsne_slider.value)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## UMAP
    """)
    return


@app.cell(hide_code=True)
def _(adata, fn, plot_out):
    umap_slider, display_umap = fn.make_dr_slider(
        adata,
        "umap",
        plot_out / "umap_unintegrated",
        ["patient", "flowcell"],
        ext="svg",
        theme_kws={"figure_size": (10, 5)},
    )
    umap_slider
    return display_umap, umap_slider


@app.cell
def _(display_umap, umap_slider):
    display_umap(umap_slider.value)
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Marker genes
    """)
    return


if __name__ == "__main__":
    app.run()
