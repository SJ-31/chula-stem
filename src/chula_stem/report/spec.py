URLS: dict = {
    "dbvar": "https://www.ncbi.nlm.nih.gov/dbvar/variants",
    "pubchem": "https://pubchem.ncbi.nlm.nih.gov/compound",
}

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
    "Database Name": 90,
    "Variant Type": 90,
    "ClinVar": 100,
}

CTABLE_RENAME: dict = {
    "Type": "CNV type",
    "Known or predicted dosage-sensitive genes": "Known/predicted dosage-sensitive genes",
    "All protein coding genes": "All genes",
    "ClinGen_report": "ClinGen",
    "source": "Database/Study records",
}

CTABLE_WIDTHS: dict = {
    "Known or predicted dosage-sensitive genes": "Known/predicted dosage-sensitive genes",
    "All protein coding genes": "All genes",
    "Type": "CNV type",
}
