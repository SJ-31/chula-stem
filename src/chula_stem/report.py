from abc import ABC, abstractmethod
from pathlib import Path
from time import sleep
from typing import override

import polars as pl
import requests
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport
from requests import Session


def ensembl_id2name(id: str) -> str | None:
    """Convert ensembl gene id to HGNC gene name

    :param id: id to convert

    :returns: name of gene according to HGNC, None if gene could not be converted
    """
    try:
        gene: pe.Gene = REL.gene_by_id(id)
        return gene.gene_name
    except ValueError:
        return None


def filter_format_vep(input: str, sep="\t"):
    wanted_cols = [
        "Gene",
        "VAF",
        "Alt_depth",
        "Loc",
        "SYMBOL",
        "Consequence",
        "CLIN_SIG",
        "Known Variants",
        "HGVSg",
    ]
    to_rename: dict = {
        "Alt_depth": "Variant Allele Depth",
        "VAF": "Variant Allele Fraction",
        "Loc": "Location",
        "CLIN_SIG": "ClinVar",
        "HGVSg": "HGVS",
        "Existing_variation": "Known Variants",
        "Consequence": "Variant Effect",
    }
    df = (
        pl.read_csv(input, separator=sep)
        .with_columns(
            pl.col("Existing_variation").map_elements(
                lambda x: x.replace("&", ", "), return_dtype=pl.String
            )
        )
        .select(wanted_cols)
        .rename(to_rename)
    )
    return df


## * Therapy db parent class


class TherapyDB(ABC):

    api_wait: float
    cache: Path
    shared_cols: tuple = ("gene", "source", "disease", "therapies")

    @abstractmethod
    def get_gene(self, gene: str, confident: bool) -> pl.DataFrame:
        pass

    @staticmethod
    def _cache_default(name: str, cache: str) -> Path:
        return (
            Path(cache)
            if Path(cache).exists()
            else Path.home().joinpath(f".cache/{name}_cache.json")
        )

    @staticmethod
    @abstractmethod
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """Retrieve the most confident and useful Evidence Item from the `get_gene` operation
        Criteria are database-specific
        :return: a df of one row
        """
        pass

    def get_genes(self, gene_list: list, confident: bool = True) -> pl.DataFrame:
        results: list = []
        previous, gene_list = self.read_cache(gene_list)
        for g in gene_list:
            current: pl.DataFrame = self.get_gene(g, confident)
            if not current.is_empty():
                results.append(current)
            sleep(self.api_wait)  # API permits 2 requests per second
        if results:
            df = pl.concat(results)
            self.write_cache(df)
            if not previous.is_empty():
                df = pl.concat([previous, df])
            return df
        if not previous.is_empty():
            return previous
        return pl.DataFrame()

    @staticmethod
    @abstractmethod
    def parse_gene_response(data: list) -> pl.DataFrame:
        """
        :param: data a list of evidence entries associated with a gene, such as the result
        of an API call
        :return: a df with at least the columns
            - gene: str
            - source: str
            - disease: list[str]
            - therapies: list[str]
        """
        pass

    def read_cache(self, gene_list: list) -> tuple[pl.DataFrame, list[str]]:
        """Retrieve entries from `gene_list` found in the cache, and filter
        `gene_list` so that only unknown entries remain
        """
        if self.cache.exists():
            previous = pl.read_json(self.cache, infer_schema_length=None)
            if not previous.is_empty():
                previous = previous.filter(pl.col("gene").is_in(gene_list))
                found: list[str] = list(
                    filter(lambda x: x not in previous["gene"], gene_list)
                )
            else:
                found = gene_list
            return previous, found
        else:
            return pl.DataFrame(), gene_list

    def write_cache(self, df: pl.DataFrame) -> None:
        df.write_json(self.cache)


## * Civic database


class Civic(TherapyDB):
    def __init__(self, cache: str = "") -> None:
        transport = AIOHTTPTransport(url="https://civicdb.org/api/graphql")
        self.api_wait: float = 0.5
        self.cache: Path = self._cache_default("civic", cache)
        self.client: Client = Client(
            transport=transport, fetch_schema_from_transport=True
        )
        self.gene_query: str = gql(
            """
        query ($entrez_symbol: String, $next_page: String) {
        gene(entrezSymbol: $entrez_symbol) {
            variants(after: $next_page) {
            pageInfo {
                hasNextPage
                endCursor
            }
            nodes {
                clinvarIds
                variantAliases
                ... on GeneVariant {
                coordinates {
                    start
                    stop
                    chromosome
                }
                }
                molecularProfiles {
                nodes {
                    name
                    evidenceItems(includeRejected: false) {
                    nodes {
                        link
                        therapies {
                        name
                        }
                        evidenceRating
                        variantOrigin
                        disease {
                            name
                        }
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
        """
        )

    @staticmethod
    def format_source(source_info: dict) -> str:
        """Get markdown-formatted url from `source_info`"""
        url = source_info.get("sourceUrl", "")
        title = source_info.get("title", "NONE")
        return f"[{title}]({url})"

    @staticmethod
    def get_loc(coordinates: dict) -> str | None:
        if not all(coordinates.values()):
            return None
        return f"{coordinates['chromosome']}:{coordinates['start']}-{coordinates['stop']}"

    @staticmethod
    @override
    def parse_gene_response(data: list) -> pl.DataFrame:
        cols: dict = {
            "loc": [],
            "molecularProfile": [],
            "therapies": [],
            "evidenceLevel": [],
            "evidenceRating": [],
            "civicLink": [],
            "origin": [],
            "disease": [],
            "source": [],
            "variantAliases": [],
            "clinvarIds": [],
        }
        for entry in data:
            loc: str = Civic.get_loc(entry["coordinates"])
            for mp in entry["molecularProfiles"]["nodes"]:
                for ev in mp["evidenceItems"]["nodes"]:
                    ther: list = [t["name"] for t in ev["therapies"]]
                    if d := ev.get("disease"):
                        disease = [d["name"]]
                    else:
                        disease = ["NA"]
                    cols["clinvarIds"].append(entry["clinvarIds"])
                    cols["variantAliases"].append(entry["variantAliases"])
                    cols["molecularProfile"].append(mp["name"])
                    cols["civicLink"].append(ev["link"])
                    cols["origin"].append(ev["variantOrigin"])
                    cols["disease"].append(disease)
                    cols["evidenceRating"].append(ev["evidenceRating"])
                    cols["evidenceLevel"].append(ev["evidenceLevel"])
                    cols["source"].append(Civic.format_source(ev["source"]))
                    cols["loc"].append(loc)
                    cols["therapies"].append(ther)

        level_mapping: dict = {None: 0, "E": 1, "D": 2, "C": 3, "B": 4, "A": 5}
        entry_df: pl.DataFrame = pl.DataFrame(cols).with_columns(
            pl.col("evidenceLevel").replace_strict(level_mapping)
        )
        return entry_df

    @staticmethod
    @override
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """
        Criteria are...
            - Entry whose associated gene variant has a known location
            - Has the highest evidence rating within the evidence set, and at least 4
            - Has the highest evidence level and at least least B (encoded as 4)
            - Has a named therapy
        """
        df = df.filter(
            (pl.col("therapies").list.len() >= 1)
            & (pl.col("evidenceRating") == pl.col("evidenceRating").max())
            & (pl.col("evidenceLevel") == pl.col("evidenceLevel").max())
            & (pl.col("loc").is_not_null())
            & (pl.col("evidenceRating") >= 4)
            & (pl.col("evidenceLevel") >= 4)
            & ((pl.col("origin") == "SOMATIC") | (pl.col("origin") == "COMBINED"))
        )
        if df.shape[0] > 1:
            return df.head(1)
        return df

    @override
    def get_gene(self, gene: str, confident: bool = True) -> pl.DataFrame:
        """Retrieve civic therapeutic data for a Gene
        By default, will only return the 'most confident' therapy, which is determined
            by the 'get_confident' function

        :param gene: gene id (HUGO format)
        :param variant_aliases: existing

        :returns:
        """
        has_next: bool = True
        query_input: dict = {"entrez_symbol": gene, "next_page": ""}
        dfs: list[pl.DataFrame] = []
        while has_next:
            response: dict = self.client.execute(
                self.gene_query, variable_values=query_input
            )
            if response["gene"] and response["gene"].get("variants"):
                variants: dict = response["gene"]["variants"]
                data: list = variants["nodes"]
                has_next = variants["pageInfo"]["hasNextPage"]
                if has_next:
                    query_input["next_page"] = variants["pageInfo"]["endCursor"]
                    dfs.append(Civic.parse_gene_response(data))
            else:
                print(f"No civic entry for gene {gene} found")
                break
        if not dfs:
            return pl.DataFrame()
        df: pl.DataFrame = pl.concat(dfs).with_columns(gene=pl.lit(gene))
        if confident:
            return Civic.get_confident(df)
        else:
            return df


## * Pandrugs 2


class PanDrugs2(TherapyDB):
    def __init__(self, cache: str = "") -> None:
        self.url: str = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/"
        self.api_wait: float = 1
        self.session: Session = Session()
        self.cache: Path = self._cache_default("pandrugs2", cache)

    @staticmethod
    @override
    def parse_gene_response(data: list) -> pl.DataFrame:
        cols: dict = {
            "status": [],
            "therapies": [],
            "therapyType": [],
            "disease": [],
            "gScore": [],
            "dScore": [],
        }
        for entry in data:
            cols["disease"].append(entry["cancer"])
            cols["dScore"].append(entry["dScore"])
            cols["status"].append(entry["status"])
            cols["gScore"].append(entry["gScore"])
            cols["therapies"].append([entry["showDrugName"]])
            cols["therapyType"].append(entry["therapy"])
        return pl.DataFrame(cols).with_columns(
            source=pl.lit("[PanDrugs2](https://pandrugs.org/#!/)")
        )

    @override
    @staticmethod
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """
        Criteria are...
            - Has the highest dScore, and dScore > 0.7
            - Has gScore > 0.6
            - Drug's status is either "APPROVED" or "CLINICAL_TRIALS"
        """
        df = df.filter(
            (pl.col("gScore") >= 0.6)
            & (pl.col("gScore") >= 0.7)
            & (
                (pl.col("status") == "APPROVED")
                | (pl.col("status") == "CLINIAL_TRIALS")
            )
        )
        df = df.filter(pl.col("dScore") == pl.col("dScore").max())

        if df.shape[0] > 1:
            get_first: list = ["source", "gene", "status"]
            make_unique: list = ["disease", "gScore", "therapyType"]
            others: list = list(
                filter(lambda x: x not in make_unique, df.columns)
            )
            others.remove("dScore")
            df = (
                df.group_by(pl.col("dScore"))
                .agg(pl.col(make_unique).flatten(), pl.col(others))
                .with_columns(
                    pl.col(get_first).list.first(),
                    pl.col(make_unique).list.unique(),
                )
            )
        return df

    @override
    def get_gene(self, gene: str, confident: bool = True) -> pl.DataFrame:
        response: requests.Response = self.session.get(
            url=self.url,
            params={
                "gene": gene,
                # The following are all required
                "directTarget": True,
                "biomarker": True,
                "pathwayMember": False,  # DO NOT permit indirect gene-drug interactions
                "geneDependency": False,
            },
        )
        if not response.ok or not (data := response.json().get("geneDrugGroup")):
            return pl.DataFrame()
        df: pl.DataFrame = PanDrugs2.parse_gene_response(data).with_columns(
            gene=pl.lit(gene)
        )
        if confident:
            return self.get_confident(df)
        else:
            return df


## * Report formatter
