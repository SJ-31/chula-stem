import marimo

__generated_with = "0.19.7"
app = marimo.App()


@app.cell
def _():
    from yte import process_yaml

    import functions as fn
    import anndata as ad

    with open("./cellranger_config.yaml", "r") as f:
        env = process_yaml(f)
        print(env)

    combined: ad.AnnData = fn.prepare_data(env)
    return (combined,)


@app.cell
def _(combined: "ad.AnnData"):
    combined
    return


if __name__ == "__main__":
    app.run()
