#!/usr/bin/env ipython

import datetime
import polars as pl
import re
from bs4 import BeautifulSoup
from bs4.element import Tag
from urllib.request import urlopen

COSMIC_URL: str = "https://cancer.sanger.ac.uk/signatures"
AET: str = "Proposed_Aetiology"


def read_page(url: str) -> str:
    with urlopen(url) as f:
        bytes = f.read()
        html = bytes.decode()
        return html


def get_artefacts(soup: BeautifulSoup, collection: str) -> pl.DataFrame:
    artefact_div: Tag = soup.find(
        lambda x: any(div.text == "Possible sequencing artefacts" for div in x.children)
    )
    artefacts = artefact_div.findChildren("a")
    if not artefacts:
        return pl.DataFrame()
    data: dict = {"Signature": [], "Link": []}
    for item in artefacts:
        data["Signature"].append(item.text)
        href = item.attrs.get("href")
        if href:
            link = f"{COSMIC_URL}/{href}"
        else:
            link = None
        data["Link"].append(link)
    df = pl.DataFrame(data).with_columns(
        pl.lit("Possible sequencing artefact").alias(AET),
        pl.lit(collection).alias("Collection"),
    )
    return df


def parse_cosmic_signature_page(source: str, url: bool = False, collection: str = ""):
    """
    Parse COSMIC signature collection page into a polars dataframe
    Last updated for page when <2024-12-25 Wed>
    """
    if source and url:
        html = read_page(source)
    else:
        with open(source, "r") as f:
            html = f.read()
    soup: BeautifulSoup = BeautifulSoup(html, "html.parser")
    signatures: list[Tag] = soup.find_all("div", {"class": "signature-card"})
    data: dict = {"Signature": [], AET: [], "Link": []}
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
            data[AET].append(pa_desription)
            if collection:
                link = f"{COSMIC_URL}/{collection.lower().replace(' ', '-')}/{name.text.lower()}"
                data["Link"].append(link)
            else:
                data["Link"].append(pl.lit(None))
    df = pl.DataFrame(data).with_columns(Collection=pl.lit(collection))
    artefacts: pl.DataFrame = get_artefacts(soup, collection)
    if not artefacts.is_empty():
        df = pl.concat([df, artefacts.select(df.columns)])
    return df


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="")
    parser.add_argument("-s", "--separator", default=",")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    urls: list = [
        "https://cancer.sanger.ac.uk/signatures/sbs/",
        "https://cancer.sanger.ac.uk/signatures/dbs/",
        "https://cancer.sanger.ac.uk/signatures/id/",
        "https://cancer.sanger.ac.uk/signatures/cn/",
        "https://cancer.sanger.ac.uk/signatures/rna-sbs/",
    ]
    print("Retrieving Cosmic signatures...")
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    default_out = f"cosmic_signatures_v3.4-{date}.csv"
    output = args["output"] if args["output"] else default_out
    all_sigs = pl.concat([parse_cosmic_signature_page(u, True) for u in urls])
    all_sigs.write_csv(output, separator=args["separator"])
