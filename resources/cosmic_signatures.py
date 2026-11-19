#!/usr/bin/env ipython

import datetime
import re
from urllib.request import urlopen

import polars as pl
from bs4 import BeautifulSoup
from bs4.element import Tag

COSMIC_URL: str = "https://cancer.sanger.ac.uk/signatures"
AET: str = "Proposed_Aetiology"


def read_page(url: str) -> str:
    with urlopen(url) as f:
        bytes = f.read()
        html = bytes.decode()
        return html


def parse_cosmic_signature_page(source: str, url: bool = False, collection: str = ""):
    """
    Parse COSMIC signature collection page into a polars dataframe
    Last updated [2026-11-19 Thu]
    """
    if source and url:
        html = read_page(source)
    else:
        with open(source, "r") as f:
            html = f.read()
    soup: BeautifulSoup = BeautifulSoup(html, "html.parser")
    signatures: list[Tag] = soup.find_all("div", {"class": "signature-card"})
    data: dict = {"Signature": [], AET: [], "Link": []}
    find_collection: list[str] = re.findall(
        r".*\| (.*) - Mutational Signatures.*", soup.title.text
    )
    if find_collection:
        collection = find_collection[0]
    for sig in signatures:
        if sig_title := sig.find("div", {"class": "signature-card-title"}):
            name: Tag | None = sig_title.find("h4")
            pa: Tag | None = sig.find("div", {"class": "signature-card-body"})
            if pa:
                data["Signature"].append(name.text)
                pa_desription = pa.text.strip().replace("Proposed Aetiology", "")
                data[AET].append(pa_desription)
                if collection:
                    link = f"{COSMIC_URL}/{collection.lower().replace(' ', '-')}/{name.text.lower()}"
                    data["Link"].append(link)
                else:
                    data["Link"].append(pl.lit(None))
    df = pl.DataFrame(data).with_columns(Collection=pl.lit(collection))
    return df


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="")
    parser.add_argument("-s", "--separator", default=",")
    args = vars(parser.parse_args())
    return args


def main() -> pl.DataFrame:
    urls: list = [
        "https://cancer.sanger.ac.uk/signatures/sbs/",
        "https://cancer.sanger.ac.uk/signatures/dbs/",
        "https://cancer.sanger.ac.uk/signatures/id/",
        "https://cancer.sanger.ac.uk/signatures/cn/",
        "https://cancer.sanger.ac.uk/signatures/rna-sbs/",
    ]
    print("Retrieving Cosmic signatures...")
    all_sigs = pl.concat(
        [parse_cosmic_signature_page(u, True) for u in urls], how="diagonal_relaxed"
    )
    return all_sigs


if __name__ == "__main__":
    args = parse_args()

    date = datetime.datetime.now().strftime("%Y-%m-%d")
    default_out = f"cosmic_signatures_v3.6-{date}.csv"
    output = args["output"] if args["output"] else default_out
    all_sigs = main()
    all_sigs.write_csv(output, separator=args["separator"])
