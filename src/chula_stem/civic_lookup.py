#!/usr/bin/env python3

from pathlib import Path

import pandera.polars as pa
import polars as pl
from attrs import define, field
from chula_stem.databases import Civic
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport


@define
class CivicVariants:
    """
    A class for performing looking up gene variants in the CIVIC database
    by their HGVSc codes
    """

    client: Client = field(
        factory=lambda: Client(
            transport=AIOHTTPTransport(url="https://civicdb.org/api/graphql"),
            fetch_schema_from_transport=True,
        )
    )
    query = gql("""
        query ($entrez_symbol: String, $next_page: String) {
        gene(entrezSymbol: $entrez_symbol) {
            variants(after: $next_page) {
            pageInfo {
                hasNextPage
                endCursor
            }
            nodes {
                hgvsDescriptions
                molecularProfiles {
                nodes {
                    name
                    evidenceItems(includeRejected: false) {
                    nodes {
                        link
                        evidenceRating
                        evidenceLevel
                        source {
                            sourceType
                            sourceUrl
                            title
                        }
                    }
                    }
                }
                }
            }
            }
        }
        }
    """)

    def _parse_variants(
        self, data: list[dict], wanted_hgvscs: list[str]
    ) -> pl.DataFrame:
        hgvscs: set = set(wanted_hgvscs)
        tmp = {
            "HGVSc": [],
            "link": [],
            "evidenceLevel": [],
            "evidenceRating": [],
            "source": [],
        }
        for vdict in data:
            if not (hgvsc := vdict.get("hgvsDescriptions")):
                continue
            intersection = set(hgvsc) & hgvscs
            if not intersection:
                continue
            hgvscs -= intersection
            for mp in vdict["molecularProfiles"]["nodes"]:
                for ev in mp["evidenceItems"]["nodes"]:
                    tmp["HGVSc"].extend(intersection)
                    tmp["in_civic"].extend([True] * len(intersection))
                    for k in tmp.keys():
                        if k == "source":
                            val = Civic.format_source(ev[k])
                        else:
                            val = ev[k]
                        tmp[k].extend([val] * len(intersection))
        if not tmp["HGVSc"]:
            return pl.DataFrame({"HGVSc": wanted_hgvscs}).with_columns(
                pl.lit(False).alias("in_civic")
            )
        df = pl.DataFrame(tmp)
        if hgvscs:
            remaining = pl.DataFrame({"HGVSc": hgvscs}).with_columns(
                pl.lit(False).alias("in_civic")
            )
            df = pl.concat([df, remaining])
        return df

    def _get_symbol_variants(self, symbol: str, hgvscs: list[str]):
        has_next: bool = True
        query_input: dict = {"entrez_symbol": symbol, "next_page": ""}
        dfs: list[pl.DataFrame] = []
        while has_next:
            response: dict = self.client.execute(
                self.query, variable_values=query_input
            )
            if response["gene"] and response["gene"].get("variants"):
                variants: dict = response["gene"]["variants"]
                data: list = variants["nodes"]
                if not data:
                    break
                df = self._parse_variants(data, hgvscs)
                hgvscs = [hg for hg in hgvscs if hg not in df["HGVSc"]]
                dfs.append(df)
                has_next = variants["pageInfo"]["hasNextPage"]
                if has_next:
                    query_input["next_page"] = variants["pageInfo"]["endCursor"]
            else:
                break
        if not dfs:
            return pl.DataFrame({"HGVSc": hgvscs}).with_columns(
                pl.lit(symbol).alias("symbol"), pl.lit(False).alias("in_civic")
            )
        df: pl.DataFrame = pl.concat(dfs, how="vertical_relaxed").with_columns(
            symbol=pl.lit(symbol)
        )
        return df

    def __call__(
        self,
        df: pl.DataFrame,
        symbol_col: str = "SYMBOL",
        hgvsc_col: str = "HGVSc",
        previous: pl.DataFrame | None = None,
    ) -> pl.DataFrame:
        schema = pa.DataFrameSchema(
            {symbol_col: pa.Column(pl.String), hgvsc_col: pa.Column(pl.String)},
            unique=[symbol_col, hgvsc_col],
        )
        df = schema.validate(df)
        if (
            previous is not None
            and "in_civic" in previous.columns
            and hgvsc_col in previous.columns
        ):
            df = df.filter(pl.col(hgvsc_col).is_in(previous[hgvsc_col]))
        results = []
        for symbol, grouped in df.group_by(symbol_col):
            assert isinstance(symbol[0], str)
            cur = self._get_symbol_variants(symbol[0], list(grouped[hgvsc_col]))
            results.append(cur)
        if previous is None:
            return pl.concat(results, how="vertical_relaxed")
        return pl.concat(results + [previous], how="vertical_relaxed")


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+")
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "-s",
        "--symbol_column",
        default="SYMBOL",
        help="Column in input tsv files containing HGNC symbols",
        action="store",
    )
    parser.add_argument(
        "-v",
        "--hgvsc_col",
        default="HGVSc",
        help="Column in input tsv files containing HGVSc string for variants",
        action="store",
    )
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
