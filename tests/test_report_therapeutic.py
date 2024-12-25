#!/usr/bin/env ipython

from collections import Counter
import re
from os import replace

import polars as pl
import polars.selectors as cs
from chula_stem.report.format import dbvar_link, get_clingen_link
from chula_stem.utils import read_facets_rds
from chula_stem.report import ReportElement, ResultsReport
from chula_stem.report.spec import URL, Rename, Widths
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, inch
from reportlab.pdfgen.canvas import Canvas

## Working on formatting the therapeutics table
table_spec = {"header_pos": (A4[0] - 5 * inch, A4[1] - inch)}
numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)

# Will have
# May have to have tables in tables to show the therapy info properly

# For small and structural variants
# TODO: will also need to do this for mutational signatures.
# CNV data is a bonus, but probably not possible
#

#
# others = df.filter(pl.col("db") != "pandrugs2")
# * Treatment
# # pubchemId
# * Evidence type
#   * Can be mutational signature or variant
# * Relevant evidence in sample (gene names or mutational signature ids)
# * Use in other cancers
# *  Study
# # Db link
#

# Is present if sig
# contributes to more than 100 mutations or exceeds 25% of the mutational activity

total_input_mutations = 55
solution_activities = "/home/shannc/Bio_SDD/chula-stem/tests/6-null-SigProfilerAssignment/Activities/Assignment_Solution_Activities.txt"
solution_activities = "/home/shannc/Bio_SDD/chula-stem/tests/6-null-SigProfilerAssignment/Activities/sol_dummy.txt"

abs_threshold: int = (
    100  # Minimum absolute number of mutations required to be considered
)
rel_threshold: float = (
    0.25  # Minimum percentage of mutational activity to be considered
)
excluded_signatures: list = []

## Signature tabl
df = pl.read_csv(solution_activities, separator="\t").select(
    ~cs.by_name(*excluded_signatures)
)
signatures: list = df.columns[1:]
sig_cols = pl.col(signatures)
samples: pl.Series = df["Samples"]
replace_expr = [pl.col(s).replace_strict({True: s, False: None}) for s in signatures]

m: pl.Series = df.select(sig_cols).sum_horizontal()
# Total number of signature mutations per sample
frequencies = df.with_columns(sig_cols / m).unpivot(
    on=signatures, index="Samples", variable_name="Signatures", value_name="Frequency"
)

sums: pl.DataFrame = pl.DataFrame({"Samples": df["Samples"], "m": m})
filtered = (
    (
        df.with_columns((sig_cols > abs_threshold) | ((sig_cols / m) > rel_threshold))
        .with_columns(replace_expr)
        .with_columns(Signatures=pl.concat_list(signatures).list.drop_nulls())
        .select(["Samples", "Signatures"])
        .join(sums, on="Samples")
        .explode("Signatures")
    )
    .join(frequencies, on=["Samples", "Signatures"])
    .group_by("Samples")
    .agg(pl.col(["Signatures", "Frequency"]), pl.col("m").first())
)

# TODO: Want to display at the top of this table the total number of sig mutations (m),
# and the fraction m/ (total number of mutations)
# Should also show the number of sigs before and after filtering
# Or maybe put those to be a separate


def make_sample_table(df) -> pl.DataFrame:
    """Format a row of the signature df (a single sample) into a df to present with
    reportlab

    The signature df has four columns: Samples, Signatures, Frequency and m
    Signatures and Frequency are list columns of the same length, containing the
    names and frequencies of signatures kept in the sample in their respective orders

    :returns:
    """
    data: dict = {"Signature": [], "Frequency": [], "Type": [], "COSMIC": []}
    return df


##
from bs4 import BeautifulSoup
from bs4.element import Tag
from urllib.request import urlopen

COSMIC_URL: str = "https://cancer.sanger.ac.uk/signatures"


def parse_cosmic_signature_page(source: str, url: bool = False, collection: str = ""):
    """
    Parse COSMIC signature collection page into a polars dataframe
    Last updated for page when <2024-12-25 Wed>
    """
    if source and url:
        with urlopen(source) as f:
            bytes = f.read()
            html = bytes.decode()
    else:
        with open(source, "r") as f:
            html = f.read()
    soup = BeautifulSoup(html, "html.parser")
    signatures: list[Tag] = soup.find_all("div", {"class": "signature-card"})
    aet: str = "Proposed_aetiology"
    data: dict = {"Signature": [], aet: [], "Link": []}
    find_collection: str = re.findall(
        r".*\| (.*) - Mutational Signatures.*", soup.title.text
    )
    if find_collection:
        collection = find_collection[0]
    for sig in signatures:
        name: Tag = sig.find("h4", {"class": "signature-card-title"})
        pa: Tag = sig.find("div", {"class": "signature-card-body"})
        if name and pa:
            data["Signature"].append(name.text)
            pa_desription = pa.text.strip().replace("Proposed Aetiology", "")
            data[aet].append(pa_desription)
            if collection:
                link = f"{COSMIC_URL}/{collection.lower().replace(' ', '-')}/{name.text.lower()}"
                data["Link"].append(link)
            else:
                data["Link"].append(pl.lit(None))
    df = pl.DataFrame(data).with_columns(Collection=pl.lit(collection))
    return df


urls: list = [
    "https://cancer.sanger.ac.uk/signatures/sbs/",
    "https://cancer.sanger.ac.uk/signatures/dbs/",
    "https://cancer.sanger.ac.uk/signatures/id/",
    "https://cancer.sanger.ac.uk/signatures/cn/",
    "https://cancer.sanger.ac.uk/signatures/rna-sbs/",
]

all_sigs = pl.concat([parse_cosmic_signature_page(u, True) for u in urls])
all_sigs.write_csv(
    "/home/shannc/Bio_SDD/chula-stem/nextflow/config/cosmicv3.4_signatures.csv"
)
