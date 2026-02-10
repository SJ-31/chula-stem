import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell
def _():
    import sys

    import anndata as ad
    import marimo as mo
    import matplotlib
    import pandas as pd
    import plotnine as gg
    from pyhere import here
    from yte import process_yaml

    workdir = here("analyses", "scrnaseq_hcc")

    sys.argv.append("test=True")
    with open(workdir / "cellranger_config.yaml", "r") as f:
        env = process_yaml(f)

    def load_file_as_module(name, location):
        from importlib import util

        spec = util.spec_from_file_location(name, location)
        module = util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    fn = load_file_as_module("functions", workdir / "functions.py")
    return ad, env, fn, mo, workdir


@app.cell
def _(ad, env, fn, workdir):
    adata = ad.read_h5ad(workdir / "test" / "filtered.h5ad")
    adata = adata[:, ~adata.var["hgnc_symbol"].isnull()]

    fn.add_saved_dr(adata, env)
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Dimensionality reduction (unintegrated)
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## T-SNE
    """)
    return


@app.cell
def _(adata, fn, workdir):
    # TODO: add this to the main notebook
    tsne_slider, display_tsne = fn.make_dr_slider(
        adata,
        "t-sne",
        workdir / "tsne_plots",
        ["patient", "type"],
        ext="svg",
        theme_kws={"figure_size": (10, 10)},
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


@app.cell
def _(adata, fn, workdir):
    # TODO: add this to the main notebook
    umap_slider, display_umap = fn.make_dr_slider(
        adata,
        "umap",
        workdir / "umap_plots",
        ["patient", "type"],
        ext="svg",
        theme_kws={"figure_size": (10, 5)},
    )
    umap_slider
    return display_umap, umap_slider


@app.cell(hide_code=True)
def _(display_umap, umap_slider):
    display_umap(umap_slider.value)
    return


if __name__ == "__main__":
    app.run()
