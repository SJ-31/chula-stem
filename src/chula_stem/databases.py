import polars as pl
from abc import ABC, abstractmethod
from pathlib import Path
from time import sleep
import requests
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport
from typing import override
from requests import Session

## * Therapy db parent class


class TherapyDB(ABC):

    api_wait: float
    name: str
    cache: Path
    shared_cols: tuple = (
        "gene",
        "source",
        "disease",
        "therapies",
        "db",
        "db_link",
    )
    all_cached: pl.DataFrame

    @abstractmethod
    def get_gene(self, gene: str, confident: bool) -> pl.DataFrame:
        pass

    @staticmethod
    def _cache_default(name: str, cache: str) -> Path:
        return (
            Path(cache) if cache else Path.home().joinpath(f".cache/{name}_cache.json")
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
        previous, gene_list = self.read_cache(gene_list)  # Get cached results
        # and filter gene list for those that weren't cached
        for g in gene_list:
            current: pl.DataFrame = self.get_gene(g, confident)
            if not current.is_empty():
                results.append(current.with_columns(db=pl.lit(self.name)))
            sleep(self.api_wait)  # API permits 2 requests per second
        if results:
            df = pl.concat(results, how="vertical_relaxed")
            if not previous.is_empty():
                df = pl.concat([previous, df], how="vertical_relaxed")
            self.write_cache(df)
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
        self.all_cached = pl.DataFrame()
        if self.cache.exists():
            previous = pl.read_json(self.cache, infer_schema_length=None)
            self.all_cached = previous
            if not previous.is_empty():
                previous = previous.filter(pl.col("gene").is_in(gene_list))
                not_found: list[str] = list(
                    filter(lambda x: x not in previous["gene"], gene_list)
                )
            else:
                not_found = gene_list
            return previous, not_found
        else:
            return pl.DataFrame(), gene_list

    def write_cache(self, df: pl.DataFrame) -> None:
        if not self.all_cached.is_empty():
            df = pl.concat([self.all_cached, df], how="vertical_relaxed")
            df = df.filter(~df.is_duplicated())
        df.write_json(self.cache)


## * Civic database

CIVIC_URL: str = "https://civicdb.org"


class Civic(TherapyDB):
    def __init__(self, cache: str = "") -> None:
        transport = AIOHTTPTransport(url="https://civicdb.org/api/graphql")
        self.name: str = "Civic"
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
                        myChemInfo {
                                pubchemCid
                            }
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
        return (
            f"{coordinates['chromosome']}:{coordinates['start']}-{coordinates['stop']}"
        )

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
        schema: dict = {
            "loc": pl.String,
            "molecularProfile": pl.String,
            "therapies": pl.List(pl.String),
            "evidenceRating": pl.Int64,
            "civicLink": pl.String,
            "origin": pl.String,
            "disease": pl.List(pl.String),
            "source": pl.String,
            "variantAliases": pl.List(pl.String),
            "clinvarIds": pl.List(pl.String),
        }
        for entry in data:
            loc: str = Civic.get_loc(entry["coordinates"])
            for mp in entry["molecularProfiles"]["nodes"]:
                for ev in mp["evidenceItems"]["nodes"]:
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
                    if ev["therapies"]:
                        ther: list = [t["name"] for t in ev["therapies"]]
                        pubchem = []
                        for t in ev["therapies"]:
                            if (
                                t
                                and (tmp := t.get("myChemInfo", {}))
                                and (id := tmp.get("pubchemCid"))
                            ):
                                pubchem.append(id)
                            else:
                                pubchem.append("NA")
                        pubchem = list(filter(lambda x: x, pubchem))
                        with_pc = [f"{x}:{y}" for x, y in zip(ther, pubchem)]
                        cols["therapies"].append(with_pc)
                    else:
                        cols["therapies"].append([])

        level_mapping: dict = {
            "": 0,
            None: 0,
            "E": 1,
            "D": 2,
            "C": 3,
            "B": 4,
            "A": 5,
        }
        entry_df: pl.DataFrame = pl.DataFrame(cols).cast(schema)
        if entry_df.shape[0] > 0:
            entry_df = entry_df.with_columns(
                pl.col("evidenceLevel").replace_strict(level_mapping),
                pl.col("civicLink").str.replace("^", CIVIC_URL).alias("db_link"),
            )
        else:
            entry_df = entry_df.cast({"evidenceLevel": pl.Int64})
        return entry_df

    @staticmethod
    @override
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """
        Criteria are...
            - Entry whose associated gene variant has a known location
            - Has the highest evidence rating within the evidence set, and at least 4
            - Has the highest evidence level and at least least B (encoded as 4)
        """
        df = df.filter(
            (pl.col("evidenceRating") == pl.col("evidenceRating").max())
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
                if not data:
                    print(f"No civic entry for gene {gene} found")
                    break
                dfs.append(Civic.parse_gene_response(data))
                has_next = variants["pageInfo"]["hasNextPage"]
                if has_next:
                    query_input["next_page"] = variants["pageInfo"]["endCursor"]
            else:
                print(f"No civic entry for gene {gene} found")
                break
        if not dfs:
            return pl.DataFrame()
        df: pl.DataFrame = pl.concat(dfs, how="vertical_relaxed").with_columns(
            gene=pl.lit(gene)
        )
        if confident:
            return Civic.get_confident(df)
        else:
            return df


## * Pandrugs 2


class PanDrugs2(TherapyDB):
    def __init__(self, cache: str = "") -> None:
        self.url: str = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/"
        self.api_wait: float = 1
        self.name: str = "PanDrugs2"
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
        schema: dict = {
            "status": pl.String,
            "therapies": pl.String,
            "therapyType": pl.String,
            "disease": pl.List(pl.String),
            "gScore": pl.Float64,
            "dScore": pl.Float64,
        }
        for entry in data:
            cols["disease"].append(entry["cancer"])
            # All diseases from pandrugs are cancers
            cols["dScore"].append(entry["dScore"])
            cols["status"].append(entry["status"])
            therapy = entry["showDrugName"]
            if entry.get("pubchemId"):
                pc = entry["pubchemId"][0]  # These are synonyms, so just take the first
                cols["therapies"].append(f"{therapy}:{pc}")
            else:
                cols["therapies"].append(f"{therapy}:NA")
            cols["gScore"].append(entry["gScore"])
            cols["therapyType"].append(entry["therapy"])
        return pl.DataFrame(cols, schema=schema).with_columns(
            db_link=pl.lit("https://pandrugs.org/#!/"),
            source=pl.lit("NA"),
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

        get_first: list = ["source", "gene", "status"]
        make_unique: list = ["disease", "gScore", "therapyType"]
        if df.shape[0] > 1:
            others: list = list(filter(lambda x: x not in make_unique, df.columns))
            others.remove("dScore")
            df = (
                df.group_by(pl.col("dScore"))
                .agg(pl.col(make_unique).flatten(), pl.col(others))
                .with_columns(
                    pl.col(get_first).list.first(),
                    pl.col(make_unique).list.unique(),
                )
            )
        elif not df.is_empty():
            to_list = list(
                filter(
                    lambda x: x not in get_first + ["dScore", "disease"],
                    df.columns,
                )
            )
            list_exprs = [pl.col(u).map_elements(lambda x: [x]) for u in to_list]
            df = df.with_columns(list_exprs)
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
            pl.lit(gene).alias("gene"),
            pl.col("therapies").map_elements(
                lambda x: [x], return_dtype=pl.List(pl.String)
            ),
        )
        if confident:
            df = self.get_confident(df)
        return df.select(sorted(df.columns))


def add_therapy_info(
    vep_out: str, civic_cache: str = "", pandrugs2_cache: str = ""
) -> pl.DataFrame:
    """Add therapeutic information

    :param vep_out: the filtered tsv from VEP reporting the effects of small and
    structural variants

    :returns: the input dataframe with additional information about therapeutic
    options for each gene.
    Specifically, adds the columns
    - "VAR_ID" to uniquely identify variants,
    - "source" database or paper source
    - "disease": disease(s) the gene is implicated in
    - "therapies": validated therapies for the gene
    """
    df: pl.DataFrame = pl.read_csv(
        vep_out, separator="\t", null_values=["NA", "."], infer_schema_length=None
    ).with_columns(
        pl.concat_str([pl.col("Loc"), pl.col("Feature")], separator="|").alias("VAR_ID")
    )
    therapy_dbs: list[TherapyDB] = [
        Civic(civic_cache),
        PanDrugs2(pandrugs2_cache),
    ]
    gene_list = list(df["SYMBOL"].unique())

    temp: list[pl.DataFrame] = []
    for db in therapy_dbs:
        # <2024-12-20 Fri>
        # BUG: setting the confidence filters off is temporary, for testing
        find_info: pl.DataFrame = db.get_genes(gene_list, False)
        if not find_info.is_empty():
            found = find_info["gene"]
            gene_list = list(filter(lambda x: x not in found, gene_list))
            temp.append(find_info.select(TherapyDB.shared_cols))
        if not gene_list:
            break
    if temp:
        all_drug_info: pl.DataFrame = pl.concat(temp, how="vertical_relaxed")
        with_therapeutics = df.join(
            all_drug_info, how="left", left_on="SYMBOL", right_on="gene"
        )
        return with_therapeutics
    return df.with_columns([pl.lit(None).alias(c) for c in TherapyDB.shared_cols])
