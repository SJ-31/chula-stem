# + Can classify variants by clinical significance
# + Metrics to include
#   * VAF (calculated b), coverage
#   * Common gene name if available
#   * Relevant therapies
#     * In this, and other cancer types
#   * VCF-level
import os

import polars as pl

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


##
## * Civic database
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport


class Civic:
    def __init__(self) -> None:
        transport = AIOHTTPTransport(url="https://civicdb.org/api/graphql")
        self.client = Client(
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
    def parse_gene_response(data: list) -> pl.DataFrame:
        cols: dict = {
            "Loc": [],
            "molecularProfile": [],
            "therapy": [],
            "evidenceRating": [],
            "source": [],
            "variantAliases": [],
            "clinvarIds": [],
        }
        for entry in data:
            loc: str = Civic.get_loc(entry["coordinates"])
            for mp in entry["molecularProfiles"]["nodes"]:
                for ev in mp["evidenceItems"]["nodes"]:
                    for ther in ev["therapies"]:
                        cols["clinvarIds"].append(entry["clinvarIds"])
                        cols["variantAliases"].append(entry["variantAliases"])
                        cols["molecularProfile"].append(mp["name"])
                        cols["evidenceRating"].append(ev["evidenceRating"])
                        cols["source"].append(Civic.format_source(ev["source"]))
                        cols["Loc"].append(loc)
                        cols["therapy"].append(ther["name"])
        entry_df: pl.DataFrame = pl.DataFrame(cols)
        return entry_df

    def get_confident(df: pl.DataFrame, existing_variation: list) -> pl.DataFrame:
        df = df.filter(
            (pl.col("evidenceRating") == pl.col("evidenceRating").max())
            & (pl.col("Loc").is_not_null())
        )
        return df

    #     return gene_df

    def get_gene(
        self, gene: str, existing_variation: list, confident: bool = True
    ) -> pl.DataFrame:
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
            has_next: bool = variants["pageInfo"]["hasNextPage"]
            if has_next:
                query_input["next_page"] = variants["pageInfo"]["endCursor"]
            dfs.append(Civic.parse_gene_response(data))
        if not dfs:
            return pl.DataFrame()
        df: pl.DataFrame = pl.concat(dfs).with_columns(Gene=pl.lit(gene))
        if confident:
            return Civic.get_confident(df, existing_variation)
        else:
            return df

    # def get_genes(self, gene_list: str)
    # TODO: this function needs to have API delays


civ = Civic()
braf = civ.get_gene("BRAF", [], False)


## * Pandrugs 2
import requests

url = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/"
allowed_status = "APPROVED"
gene = "KDR"
response = requests.get(
    url=url,
    params={
        "gene": gene,
        # The following are all required
        "directTarget": True,
        "cancerDrugStatus": allowed_status,
        "biomarker": True,
        "pathwayMember": False,
        "geneDependency": False,
    },
)
data = response.json()
data["geneDrugGroup"][0]
