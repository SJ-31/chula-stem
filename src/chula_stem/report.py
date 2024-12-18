from abc import ABC, abstractmethod
from pathlib import Path
from time import sleep
from typing import override

import polars as pl
import polars.selectors as cs
import requests
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, SimpleDocTemplate
from requests import Session

from chula_stem.callset_qc import IMPACT_MAP


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
    all_cached: pl.DataFrame

    @abstractmethod
    def get_gene(self, gene: str, confident: bool) -> pl.DataFrame:
        pass

    @staticmethod
    def _cache_default(name: str, cache: str) -> Path:
        return (
            Path(cache)
            if cache
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
        previous, gene_list = self.read_cache(gene_list)  # Get cached results
        # and filter gene list for those that weren't cached
        for g in gene_list:
            current: pl.DataFrame = self.get_gene(g, confident)
            if not current.is_empty():
                results.append(current)
            sleep(self.api_wait)  # API permits 2 requests per second
        if results:
            df = pl.concat(results)
            if not previous.is_empty():
                df = pl.concat([previous, df])
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
            df = pl.concat([self.all_cached, df])
            df = df.filter(~df.is_duplicated())
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
            cols["therapies"].append(entry["showDrugName"])
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

        get_first: list = ["source", "gene", "status"]
        make_unique: list = ["disease", "gScore", "therapyType"]
        if df.shape[0] > 1:
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
            gene=pl.lit(gene)
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
    options for each gene
    """
    df: pl.DataFrame = pl.read_csv(
        vep_out, separator="\t", null_values=["NA", "."], infer_schema_length=None
    ).with_columns(
        pl.concat_str([pl.col("Loc"), pl.col("Feature")], separator="|").alias(
            "VAR_ID"
        )
    )
    therapy_dbs: list[TherapyDB] = [
        Civic(civic_cache),
        PanDrugs2(pandrugs2_cache),
    ]
    gene_list = list(df["SYMBOL"].unique())

    temp: list[pl.DataFrame] = []
    for db in therapy_dbs:
        find_info: pl.DataFrame = db.get_genes(gene_list, True)
        if not find_info.is_empty():
            found = find_info["gene"]
            gene_list = list(filter(lambda x: x not in found, gene_list))
            temp.append(find_info.select(TherapyDB.shared_cols))
        if not gene_list:
            break
    if temp:
        all_drug_info: pl.DataFrame = pl.concat(temp)
        with_therapeutics = df.join(
            all_drug_info, how="left", left_on="SYMBOL", right_on="gene"
        )
        return with_therapeutics
    return df.with_columns([pl.lit(None).alias(c) for c in TherapyDB.shared_cols])


## * Report formatter
def style_cells(
    start: tuple,
    ncols: int = 0,
    nrows: int = 0,
    fontname: str = "",
    fontsize: str = "",
    align: str = "",
    background: str = "",
    valign: str = "",
    textcolor=None,
    underline: tuple = (),
) -> list:
    """
    Coordinates for table style are given as (column, row)
    """
    parameter_map: dict = {
        fontname: "FONTNAME",
        fontsize: "FONTSIZE",
        textcolor: "TEXTCOLOR",
        background: "BACKGROUND",
        align: "ALIGN",
        valign: "VALIGN",
        underline: "LINEBELOW",
    }
    if ncols and nrows:
        end: tuple = start[0] + ncols - 1, start[1] + nrows - 1
    else:
        end = (-1, -1)

    def style_helper(format: str, value) -> tuple:
        if isinstance(value, tuple):
            return (format, start, end) + value
        else:
            return (format, start, end, value)

    styles: list = []
    for param, name in parameter_map.items():
        if param:
            styles.append(style_helper(name, param))
    return styles


def add_pstyles(data: pl.DataFrame | list, style: ParagraphStyle | dict) -> list:
    """Helepr for adding paragraph styles to lists or data in dfs

    :param style: A single style which is then applied to all data.
    Alternatively, a map of column_index -> style specifying styles to
    apply to specific columns. A key for 'None' is the default and applied to columns
    not explicitly given

    :returns:
    """
    if isinstance(style, dict):
        if not style.get(None):
            raise ValueError("A default key `None` must be provided!")

        style_fn = lambda x, index=None: Paragraph(str(x), style.get(index))
    else:
        style_fn = lambda x, index=None: Paragraph(str(x), style)

    if isinstance(data, pl.DataFrame):
        return [
            [Paragraph(str(s), style.get(i)) for i, s in enumerate(row)]
            for row in data.iter_rows()
        ]
    else:
        return [style_fn(row, i) for i, row in enumerate(data)]


## Report configuration
VTABLE_RENAME: dict = {
    "SYMBOL": "Gene",
    "VAF": "Variant Allele Frequency",
    "Alt_depth": "Variant Read Support",
    "Loc": "Locus",
    "HGVSc": "HGVS",
    "Existing_variation": "Database Name",
    "Consequence": "Variant Type",
    "CLIN_SIG": "ClinVar",
}
VTABLE_COL_WIDTHS: dict = {
    "Gene": 60,
    "Variant Allele Frequency": 55,
    "Variant Read Support": 45,
    "Locus": 80,
    "HGVS": 90,
    "Variant Type": 90,
    "ClinVar": 100,
}

## ** Report class


class ResultsReport:
    def __init__(
        self,
        filename: str,
        pandrugs2_cache: str,
        civic_cache: str,
        vep_small: str,
        vep_sv: str,
    ) -> None:
        self.doc: SimpleDocTemplate = SimpleDocTemplate(
            filename, rightMargin=2, leftMargin=2
        )
        self.civic_cache = civic_cache
        self.pandrugs2_cache = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}

        self._format_vep(vep_small, "small")
        self._format_vep(vep_sv, "sv")

    def _format_vep(self, df_path: str, variant_class: str) -> None:
        """Format and filter vep output into a dataframe with values ready to write
        into a reportlab table
        """
        var_col: str = "Database Name"  # New column for 'Existing_variation'
        with_therapeutics: pl.DataFrame = add_therapy_info(
            df_path, self.civic_cache, self.pandrugs2_cache
        )
        wanted_cols: list = list(VTABLE_RENAME.values())
        with_therapeutics = (
            with_therapeutics.drop("Gene")
            .rename(VTABLE_RENAME)
            .with_columns(
                pl.col("HGVS").str.extract(r".*:(.*)$", 1),
                pl.col("Variant Type").str.replace("_variant$", ""),
                (pl.col("Variant Allele Frequency").cast(pl.Float64) * 100)
                .round(3)
                .map_elements(lambda x: f"{x}%", return_dtype=pl.String),
            )
            .with_columns(
                cs.by_dtype(pl.String)
                .fill_null("-")
                .str.replace_all("&", ",\n", literal=True)
                .str.replace_all("_", " ", literal=True),
            )
        )
        confident = (
            with_therapeutics.filter(
                (pl.col(var_col).is_not_null()) & (pl.col("SOMATIC") == "1")
            )
            .with_columns(
                impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP)
            )
            .sort(pl.col("impact_score"), descending=True)
        )
        others = confident.filter(~pl.col("VAR_ID").is_in(confident["VAR_ID"]))
        self.data["relevant"][variant_class] = confident.select(wanted_cols)
        self.data["nonrelevant"][variant_class] = others.select(wanted_cols)
        self.data["all"][variant_class] = with_therapeutics.select(wanted_cols)
