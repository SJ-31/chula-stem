# + Can classify variants by clinical significance
# + Metrics to include
#   * VAF (calculated b), coverage
#   * Common gene name if available
#   * Relevant therapies
#     * In this, and other cancer types
#   * VCF-level
from abc import ABC, abstractmethod
from pathlib import Path
from time import sleep
from typing import override

import polars as pl
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport

# import pyensembl as pe

# os.environ["PYENSEMBL_CACHE_DIR"] = "/home/shannc/Bio_SDD/.cache"
# REL: pe.EnsemblRelease = pe.EnsemblRelease()

# def filter_classify_cnv()


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

    @abstractmethod
    def get_gene(self, gene: str, confident: bool) -> pl.DataFrame:
        pass

    def get_genes(self, gene_list: str, confident: bool = True) -> pl.DataFrame:
        results: list = []
        results[0], gene_list = self.read_cache(gene_list)
        for g in gene_list:
            results.append(self.get_gene(g, confident))
            sleep(self.api_wait)  # API permits 2 requests per second
        df = pl.concat(results)
        self.write_cache(df)
        return df

    @staticmethod
    @abstractmethod
    def parse_gene_response(data: list) -> pl.DataFrame:
        pass

    def read_cache(self, gene_list) -> tuple[pl.DataFrame, list[str]]:
        """Retrieve entries from `gene_list` found in the cache, and filter
        `gene_list` so that only unknown entries remain
        """
        if self.cache.exists():
            previous: pl.DataFrame = pl.read_csv(
                self.cache, separator="\t", null_values="NA"
            ).filter(pl.col("Gene").is_in(gene_list))
            found: list[str] = list(
                filter(lambda x: x in previous["Gene"], gene_list)
            )
            return previous, found
        else:
            return pl.DataFrame(), gene_list

    def write_cache(self, df: pl.DataFrame) -> None:
        if not self.cache.exists():
            df.write_csv(self.cache, separator="\t", null_value="NA")
        else:
            previous: pl.DataFrame = pl.read_csv(
                self.cache, separator="\t", null_values="NA"
            )
            write: pl.DataFrame = pl.concat([previous, df]).unique("Gene")
            write.write_csv(self.cache, separator="\t", null_value="NA")


## * Civic database


class Civic(TherapyDB):
    def __init__(self, cache: str = "") -> None:
        transport = AIOHTTPTransport(url="https://civicdb.org/api/graphql")
        self.api_wait: float = 0.5
        self.cache: Path = (
            Path(cache)
            if cache
            else Path.home().joinpath(".cache/civic_genes.tsv")
        )
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
            "Loc": [],
            "molecularProfile": [],
            "therapies": [],
            "evidenceLevel": [],
            "evidenceRating": [],
            "civicLink": [],
            "source": [],
            "variantAliases": [],
            "clinvarIds": [],
        }
        for entry in data:
            loc: str = Civic.get_loc(entry["coordinates"])
            for mp in entry["molecularProfiles"]["nodes"]:
                for ev in mp["evidenceItems"]["nodes"]:
                    ther: list = [t["name"] for t in ev["therapies"]]
                    cols["clinvarIds"].append(entry["clinvarIds"])
                    cols["variantAliases"].append(entry["variantAliases"])
                    cols["molecularProfile"].append(mp["name"])
                    cols["civicLink"].append(ev["link"])
                    cols["evidenceRating"].append(ev["evidenceRating"])
                    cols["evidenceLevel"].append(ev["evidenceLevel"])
                    cols["source"].append(Civic.format_source(ev["source"]))
                    cols["Loc"].append(loc)
                    cols["therapies"].append(ther)

        level_mapping: dict = {"E": 1, "D": 2, "C": 3, "B": 4, "A": 5}
        entry_df: pl.DataFrame = pl.DataFrame(cols).with_columns(
            pl.col("evidenceLevel").replace_strict(level_mapping)
        )
        return entry_df

    @staticmethod
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """Retrieve the most confident and useful Evidence Item from the `get_gene` operation
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
            & (pl.col("Loc").is_not_null())
            & (pl.col("evidenceRating") >= 4)
            & (pl.col("evidenceLevel") >= 4)
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
            if not response["gene"]:
                break
            variants: dict = response["gene"]["variants"]
            data: list = variants["nodes"]
            has_next = variants["pageInfo"]["hasNextPage"]
            if has_next:
                query_input["next_page"] = variants["pageInfo"]["endCursor"]
            dfs.append(Civic.parse_gene_response(data))
        if not dfs:
            return pl.DataFrame()
        df: pl.DataFrame = pl.concat(dfs).with_columns(Gene=pl.lit(gene))
        if confident:
            return Civic.get_confident(df)
        else:
            return df


## * Pandrugs 2
import requests

url = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/"
allowed_status = ["APPROVED", "CLINICAL_TRIALS"]
gene = "KDR"
response = requests.get(
    url=url,
    params={
        "gene": gene,
        # The following are all required
        "directTarget": True,
        "cancerDrugStatus": allowed_status,
        "biomarker": True,
        "pathwayMember": False,  # DO NOT permit indirect gene-drug interactions
        "geneDependency": False,
    },  # DO permit gene-drug interactions where the gene is a gene dependency
)


class Pandrugs(TherapyDB):
    def __init__(self) -> None:
        pass

    @staticmethod
    @override
    def parse_gene_response(data: list) -> pl.DataFrame:
        cols: dict = {"status": [], "therapies": [], "source": []}
        return

    @override
    def get_gene(self, gene: str, confident: bool = True) -> pl.DataFrame:
        pass
        # TODO: get the Gene here pl.col("Gene") = pl.lit(gene)


data = response.json()
data["geneDrugGroup"][0]
