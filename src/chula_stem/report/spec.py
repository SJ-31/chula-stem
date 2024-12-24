from dataclasses import dataclass


@dataclass
class URL:
    dbvar = "https://www.ncbi.nlm.nih.gov/dbvar/variants"
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/compound"


@dataclass
class Widths:
    sv_snp = {
        "Gene": 60,
        "Variant Allele Frequency": 55,
        "Variant Read Support": 45,
        "Locus": 80,
        "HGVS": 90,
        "Database Name": 90,
        "Variant Type": 90,
        "ClinVar": 100,
    }
    cnv = {
        "Known or predicted dosage-sensitive genes": "Known/predicted dosage-sensitive genes",
        "All protein coding genes": "All genes",
        "Type": "CNV type",
    }
    therapy = {
        "Therapy": 50,
        "PubChemId": 30,
        "Evidence category": 70,
        "Relevant evidence": 90,
        "Relevant cancers": 100,
        "Study": 50,
        "Database source": 40,
    }


@dataclass
class Rename:
    sv_snp = {
        "SYMBOL": "Gene",
        "VAF": "Variant Allele Frequency",
        "Alt_depth": "Variant Read Support",
        "Loc": "Locus",
        "HGVSc": "HGVS",
        "Existing_variation": "Database Name",
        "Consequence": "Variant Type",
        "CLIN_SIG": "ClinVar",
    }
    cnv = {
        "Type": "CNV type",
        "Known or predicted dosage-sensitive genes": "Known/predicted dosage-sensitive genes",
        "All protein coding genes": "All genes",
        "ClinGen_report": "ClinGen",
        "source": "Database/Study records",
    }
    therapy = {
        "Gene": "Relevant evidence",
        "disease": "Relevant cancers",
        "source": "Study source",
        "db_link": "Database source",
    }
