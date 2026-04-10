#!/usr/bin/env ipython

from pathlib import Path

import dash_cytoscape as cyto
import pandera.polars as pa
import polars as pl
from dash import Dash, Input, Output, callback, html
from pyhere import here
from yte import process_yaml

workdir = here("analyses", "scrnaseq_hcc")


with open(workdir / "cellranger_config.yaml", "r") as f:
    env = process_yaml(f)

de_graphs = Path(env["outdir"]) / "annotations" / "de_plots"

stylesheet = [
    {"selector": ".highlight", "style": {"background-color": "green"}},
    {"selector": ".other", "style": {"background-color": "gray"}},
    {"selector": "node", "style": {"content": "data(label)"}},
]


def get_df(file) -> pl.DataFrame:
    return (
        pl.read_csv(file)
        .with_row_index()
        .rename({"index": "id", "name": "label"})
        .with_columns(
            pl.when(pl.col("is_de"))
            .then(pl.lit("highlight"))
            .otherwise(pl.lit("other"))
            .alias("classes")
        )
    )


def get_edges(file) -> pl.DataFrame:
    return (
        pl.read_csv(file)
        .rename({"from": "source", "to": "target", "elabel": "label"})
        .unique(["source", "target"])
    )


test_nodes = pl.concat(
    [
        get_df(de_graphs / "leiden_res0.01_0_vs_Rest_fi_clusters_nodes.csv"),
        get_df(de_graphs / "leiden_res0.01_1_vs_Rest_fi_clusters_nodes.csv"),
    ]
).unique("id")

test_edges = pl.concat(
    [
        get_edges(de_graphs / "leiden_res0.01_0_vs_Rest_fi_clusters_edges.csv"),
        get_edges(de_graphs / "leiden_res0.01_1_vs_Rest_fi_clusters_edges.csv"),
    ]
).unique(["source", "target"])


def get_cyto_elements(node_df: pl.DataFrame, edge_df: pl.DataFrame) -> list[dict]:
    shared = {
        "label": pa.Column(pl.String, required=False),
        "classes": pa.Column(pl.String, required=False),
        "grabbable": pa.Column(pl.Boolean, required=False),
        "selectable": pa.Column(pl.Boolean, required=False),
        "locked": pa.Column(pl.Boolean, required=False),
    }
    schema = pa.DataFrameSchema(
        {"id": pa.Column(pl.String, unique=True)} | shared, coerce=True
    )
    node_df = node_df.select([c for c in node_df.columns if c in schema.columns])
    node_df = schema.validate(node_df)
    edge_schema = pa.DataFrameSchema(
        {"source": pa.Column(pl.String), "target": pa.Column(pl.String)} | shared,
        unique=["source", "target"],
        coerce=True,
    )
    edge_df = edge_df.select([c for c in edge_df.columns if c in edge_schema.columns])
    edge_df = edge_schema.validate(edge_df)
    results: list[list[dict]] = []
    node_df = node_df.filter(
        pl.col("id").is_in(edge_df["source"]) | pl.col("id").is_in(edge_df["target"])
    )
    for df, t in zip((node_df, edge_df), ("nodes", "edges")):
        cur = []
        for row in df.iter_rows(named=True):
            if t == "nodes":
                element = {"data": {"id": row["id"]}}
            else:
                element = {"data": {"source": row["source"], "target": row["target"]}}

            for k, v in row.items():
                if v is None or k in {"id", "source", "target"}:
                    continue
                if k in ["grabbable", "selectable", "locked", "classes"]:
                    element[k] = v
                else:
                    element["data"][k] = v

            cur.append(element)
        results.append(cur)

    return results[0] + results[1]


els = get_cyto_elements(test_nodes, test_edges)

# TODO: add in callbacks to display more information about each node and edge
# TODO: highlight the nearest neighbors, masking out other nodes and edges
# TODO: wrap this up in a standalone app

app = Dash()
app.layout = html.Div(
    [
        cyto.Cytoscape(
            id="graph",
            elements=els,
            style={"width": "100%", "height": "800px"},
            stylesheet=stylesheet,
            layout={
                "name": "cose",
                "nodeOverlap": 10,
                "gravity": 3,
                "nodeRepulsion": 100,
                "componentSpacing": 100,
            },
        ),
        html.Pre(id="node_data"),
    ],
)

# TODO: see https://dash.plotly.com/cytoscape/callbacks for how to do the styling
# pretty straightforward. Only difficulty is in determining the new edges


# TODO: figure out what the second arg to callback can take
# It's the property of the output element.
# `children` is commonly used because it's the property of HTML components to display new text
@callback(Output("node_data", "children"), Input("graph", "tapNodeData"))
def display_node_data(data):
    return "Gene: " + data["label"]


if __name__ == "__main__":
    app.run(debug=True, port=8021)
