import marimo

__generated_with = "0.17.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import polars as pl
    import scirpy as ir

    import fns as fn
    from fns import CONFIG, get_airrs
    return CONFIG, fn, get_airrs, ir, mo, np, pl


@app.cell
def _(CONFIG):
    CONFIG
    return


@app.cell
def _(get_airrs):
    pdac82 = get_airrs()
    return (pdac82,)


@app.cell
def _(pdac82):
    ori = pdac82["ori"]
    mixcr = pdac82["mixcr"]
    print(ori.obsm["airr"].type.show())
    return mixcr, ori


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# PDAC82 comparison""")
    return


@app.cell
def _(ir, ori, pl):
    ori_cdr3 = ir.get.airr(ori, ["complete_vdj", "cdr3"])
    pl.from_pandas(ori_cdr3)
    return


@app.cell
def _(ir, mixcr, pl):
    pl.from_pandas(ir.get.airr(mixcr, ["complete_vdj", "cdr3"]))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Comparison of null counts """)
    return


@app.cell(hide_code=True)
def _(fn, mixcr, ori):
    fn.compare_null_counts(
        [ori, mixcr],
        [
            "cdr3",
            "complete_vdj",
            "cdr1",
            "cdr2",
            "sequence",
            "v_call",
            "j_call",
            "d_call",
            "junction",
        ],
        ["bd", "mixcr"],
    )
    return


@app.cell(hide_code=True)
def _(mixcr, np):
    mchains = mixcr.obs["chain_pairing"].value_counts()
    print(mchains)
    print(
        mchains["single pair"]
        + np.sum(
            np.array(
                [
                    mchains[v] * 2
                    for v in ["extra VJ", "extra VDJ", "multichain", "two full chains"]
                ]
            )
        )
        + mchains["ambiguous"]
    )
    return


if __name__ == "__main__":
    app.run()
