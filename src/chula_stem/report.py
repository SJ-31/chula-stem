# + Can classify variants by clinical significance
# + Metrics to include
#   * VAF (calculated b), coverage
#   * Common gene name if available
#   * Relevant therapies
#     * In this, and other cancer types
#   * VCF-level
import os

import polars as pl
import pyensembl as pe

os.environ["PYENSEMBL_CACHE_DIR"] = "/home/shannc/Bio_SDD/.cache"
REL: pe.EnsemblRelease = pe.EnsemblRelease()


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
