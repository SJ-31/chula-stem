from __future__ import annotations

import copy
import importlib.resources as res
import os
from abc import ABC, abstractmethod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from time import sleep
from typing import override

import click
import duckdb
import pandera.polars as pa
import polars as pl
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport

# * Utilities


@click.command()
@click.option(
    "-i",
    "--input-files",
    required=True,
    help="Path to directory containing input files",
)
@click.option(
    "-c",
    "--gene-col",
    default="Gene",
    show_default=True,
    help="For dataframe input, the column containing gene names",
)
@click.option(
    "-f",
    "--filter-confident",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether or not to filter confident therapy entries",
)
@click.option(
    "-n",
    "--not-in-db",
    default="not_in_databases.txt",
    show_default=True,
    help="File to write genes that weren't found",
)
@click.option(
    "-i",
    "--ignore_failed_file",
    required=True,
    help="If supplied, will ignore all genes in this file",
)
@click.option(
    "-e", "--cache", required=True, help="DB file to cache data", default=None
)
@click.option(
    "-d",
    "--databases",
    required=True,
    help="Comma-delimited string specifying which databases to include in the lookup",
    default=None,
    type=str,
)
def therapy_cache_from_files(
    input_files: str,
    not_in_db: str,
    filter_confident: bool,
    cache: str,
    gene_col: str = "Genes",
):
    genes: set = set()
    for path in Path(input_files).rglob("*"):
        path: Path
        if path.is_dir():
            continue
        if path.suffix == ".tsv":
            df = pl.read_csv(path, separator="\t", infer_schema_length=None)
            genes |= set(df[gene_col])
        elif path.suffix == ".csv":
            df = pl.read_csv(path, infer_schema_length=None)
            genes |= set(df[gene_col])
        else:
            with open(path, "r") as f:
                genes |= set(f.read().splitlines())
    _ = get_therapy_df(
        genes,
        cache=Path(cache),
        not_in_db_file=not_in_db,
        filter_confident=filter_confident,
        lookup=True,
    )


def get_therapy_df(
    genes: set,
    cache: Path,
    dbs: list[str] = ["civic"],
    not_in_db_file: str = "",
    filter_confident=False,
    attempt_failed=True,
    lookup: bool = False,
) -> pl.DataFrame:
    """Get therapeutic information

    :param genes: set of genes to look up

    :returns: the input dataframe with additional information about therapeutic
    options for each gene.
    The result is a dataframe with the following columns
    - "gene"
    - "source" database or literature source
    - "disease": disease(s) the gene is implicated in
    - "therapies": validated therapies for the gene
    - "db": Database name
    - "db_link": link to page of association in database
    """
    db_mapping = {"civic": Civic, "open_targets": OpenTargets}
    therapy_dbs: list[TherapyDB] = [db_mapping[db](cache) for db in dbs]
    temp: list[pl.DataFrame] = []
    genes_copy: set = copy.copy(genes)
    if not attempt_failed and not_in_db_file and os.path.exists(not_in_db_file):
        with open(not_in_db_file, "r") as f:
            failed_previously: list = f.read().splitlines()
            genes = genes - set(failed_previously)
    for db in therapy_dbs:
        find_info: pl.DataFrame = db.get_genes(genes, filter_confident, lookup=lookup)
        TherapyDB.schema_shared_pa.validate(find_info)
        if not find_info.is_empty():
            genes -= set(find_info["gene"])
            temp.append(find_info)
        if not genes:
            break
    if temp:
        all_found: pl.DataFrame = pl.concat(temp, how="vertical_relaxed")
        not_in_db: set = genes_copy - set(all_found["gene"])
        return all_found
    else:
        not_in_db = genes_copy
    if not_in_db_file:
        with open(not_in_db_file, "a") as f:
            f.write("\n".join(not_in_db))
    return pl.DataFrame(schema=TherapyDB.schema_shared())


def no_gene_df(gene: str) -> pl.DataFrame:
    """Dataframe to return when `gene` isn't found"""
    return pl.DataFrame(
        {
            "gene": [gene],
            "found": [False],
        }
    )


# * Therapy db parent class


class TherapyDB(ABC):
    api_wait: float
    name: str
    schema: dict[str, str]
    shared_cols: tuple = (
        "gene",
        "source",
        "disease",
        "therapies",
        "db",
        "db_link",
    )
    schema_shared_pa: pa.DataFrameSchema = pa.DataFrameSchema(
        {
            "gene": pa.Column(str),
            "source": pa.Column(str),
            "disease": pa.Column(pl.List(pl.String)),
            "therapies": pa.Column(pl.List(pl.String)),
            "db": pa.Column(str),
            "db_link": pa.Column(str),
        }
    )
    batched: bool
    db: duckdb.DuckDBPyConnection

    @staticmethod
    def schema_shared() -> dict:
        return {k: v.dtype.type for k, v in TherapyDB.schema_shared_pa.columns.items()}

    def get_gene(self, gene: str, confident: bool) -> pl.DataFrame:
        raise NotImplementedError()

    @staticmethod
    def symbol2ensembl(genes: set[str]) -> pl.DataFrame:
        mapping: pl.DataFrame = pl.read_csv(
            res.files("data") / "ensembl_gene_data_116.csv"
        )
        return mapping.filter(
            (pl.col("HGNC symbol").is_in(genes))
            & (pl.col("Gene Stable ID").is_not_null())
        ).rename({"HGNC symbol": "gene", "Gene Stable ID": "ensembl"})

    @staticmethod
    def get_confident(df: pl.DataFrame) -> pl.DataFrame:
        """Retrieve the most confident and useful Evidence Item from the `get_gene` operation
        Criteria are database-specific
        :return: a df of one row
        """
        raise NotImplementedError()

    def get_genes(
        self, genes: set, confident: bool = True, lookup: bool = False
    ) -> pl.DataFrame:
        schema_concat = ",\n    ".join(
            f"{name} {dtype}" for name, dtype in self.schema.items()
        )
        self.db.execute(
            f"CREATE TABLE IF NOT EXISTS {self.name} (\n    {schema_concat}\n)"
        )
        previous: pl.DataFrame = self.db.query(f"SELECT * FROM {self.name}").pl()
        if not lookup:
            return previous.filter(pl.col("gene").is_in(genes))
        if not previous.is_empty():
            genes = genes - set(previous["gene"])
        if not self.batched:
            results: list = []
            for g in genes:
                try:
                    current: pl.DataFrame = self.get_gene(g, confident)
                    if not current.is_empty():
                        results.append(current.with_columns(db=pl.lit(self.name)))
                    sleep(self.api_wait)  # API permits 2 requests per second
                except TimeoutError:
                    print(f"WARNING: received timeout error for gene {g}, skipping...")
            if results:
                df = pl.concat(results, how="vertical_relaxed")
            else:
                df = pl.DataFrame(schema=self.schema_shared())
        else:
            df = self.get_genes_batched(genes, confident=confident)
        if not df.is_empty():
            df = df.select(self.schema.keys())
            self.db.execute(f"INSERT INTO {self.name} SELECT * FROM df")
        if not previous.is_empty():
            df = pl.concat([previous, df], how="vertical_relaxed")
        df = df.select(self.schema_shared_pa.columns.keys())
        return df

    def get_genes_batched(self, genes: set, confident: bool = True) -> pl.DataFrame:
        raise NotImplementedError()

    @staticmethod
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
        raise NotImplementedError()


# * Open targets


# * Civic database

CIVIC_URL: str = "https://civicdb.org"


class Civic(TherapyDB):
    batched = False
    name: str = "Civic"
    api_wait = 0.5
    schema: dict[str, str] = {
        "loc": "VARCHAR",
        "molecularProfile": "VARCHAR",
        "therapies": "VARCHAR[]",
        "evidenceLevel": "BIGINT",
        "evidenceRating": "BIGINT",
        "civicLink": "VARCHAR",
        "origin": "VARCHAR",
        "disease": "VARCHAR[]",
        "source": "VARCHAR",
        "variantAliases": "VARCHAR[]",
        "clinvarIds": "VARCHAR[]",
        "db_link": "VARCHAR",
        "gene": "VARCHAR",
        "db": "VARCHAR",
    }

    def __init__(self, cache: str | Path = "") -> None:
        transport = AIOHTTPTransport(url="https://civicdb.org/api/graphql")
        self.client: Client = Client(
            transport=transport, fetch_schema_from_transport=True
        )
        self.db: duckdb.DuckDBPyConnection = duckdb.connect(cache)

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
    def get_loc(coordinates: dict) -> str:
        if not all(coordinates.values()):
            return "NA"
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
        return df


class OpenTargets(TherapyDB):
    batched = True
    schema: dict[str, str] = {
        "therapies": "VARCHAR[]",
        "disease": "VARCHAR[]",
        "db_link": "VARCHAR",
        "gene": "VARCHAR",
        "db": "VARCHAR",
    }

    def __init__(self, cache, data_dir: Path | str | None = None) -> None:
        """
        Parameters
        data_dir : Path | str | None
            Optional path to directory containing pre-downloaded
            Open Targets data.
            Must have two directories with parquet files:
            `clinical_target` and `clinical_report`
        """
        self.db: duckdb.DuckDBPyConnection = duckdb.connect(cache)
        self.data_dir: Path | None = (
            Path(data_dir) if isinstance(data_dir, str) else data_dir
        )

    def get_data(self) -> tuple[duckdb.DuckDBPyRelation, duckdb.DuckDBPyRelation]:
        """
        Return PyRelation objects of Open Targets' clinical report and
        clinical targets datasets, downloading them to a temporary
        directory if necessary

        https://platform.opentargets.org/downloads
        """
        if self.data_dir is None:
            tmpdir = TemporaryDirectory()
            data_dir = Path(tmpdir.name)
            cmd_args = [
                "wget",
                "--recursive",
                "--no-parent",
                "--no-host-directories",
                "--cut-dirs",
                "6",
            ]
            _ = run(
                cmd_args
                + [
                    "ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/clinical_report",
                    "-P",
                    tmpdir.name,
                ],
                shell=True,
                check=True,
            )
            _ = run(
                cmd_args
                + [
                    "ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/clinical_target",
                    "-P",
                    tmpdir.name,
                ],
                shell=True,
                check=True,
            )
        else:
            data_dir = self.data_dir
        report = duckdb.read_parquet(f"{data_dir}/clinical_report/*.parquet")
        target = duckdb.read_parquet(f"{data_dir}/clinical_target/*.parquet")
        return report, target

    def get_genes_batched(self, genes: set, confident: bool = True) -> pl.DataFrame:
        report, target = self.get_data()
        keep = self.symbol2ensembl(genes)
        target_filtered: pl.DataFrame = target.filter(
            "targetId IN (SELECT ensembl FROM keep)"
        ).pl()
        target_filtered = (
            target_filtered.filter(
                pl.col("diseases")
                .list.filter(
                    pl.element()
                    .struct.field("diseaseFromSource")
                    .str.contains("tumor|neoplasm|cancer")
                )
                .list.len()
                > 0
            )
            .join(keep, left_on="targetId", right_on="ensembl")
            .with_columns(
                db_link="https://platform.opentargets.org/target/" + pl.col("targetId"),
                db="Open Targets",
            )
        )
        if target_filtered.is_empty():
            return pl.DataFrame()
        target_filtered = target_filtered.with_columns("" + pl.col(""))
        unique_reports = target_filtered["clinicalReportIds"].explode().value_counts()
        report = report.filter("id IN (SELECT clinicalReportIds FROM unique_reports)")
        if confident:
            report = report.filter(
                "clinicalStage IN ['APPROVAL', 'PHASE_4', 'PREAPPROVAL', 'PHASE_3']"
            )
        report_filtered: pl.DataFrame = report.pl().rename({"id": "report_id"})
        if report_filtered.is_empty():
            return (
                target_filtered.with_columns(pl.lit(None).alias("source"))
                .group_by("gene")
                .agg(pl.all())
            )
        target_filtered = (
            target_filtered.explode("clinicalReportIds")
            .join(report_filtered, left_on="clinicalReportIds", right_on="report_id")
            .select()
        )
        # TODO: make sure that the schema is the same as civic
        # and write at the top of therapy db!!!
        return


def add_therapy_info(
    vep_out: str = "",
    civic_cache: str = "",
    pandrugs2_cache: str = "",
    filter_confident=False,
    gene_list: list | None = None,
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
        find_info: pl.DataFrame = db.get_genes(gene_list, filter_confident)
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
