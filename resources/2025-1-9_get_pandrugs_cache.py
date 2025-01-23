#!/usr/bin/env ipython
import polars as pl

# <2025-01-09 Thu> Look up genes from PDAC with pandrugs2 database
# Have to do this because of some error with the pandrugs2 lookups
from chula_stem.databases import PanDrugs2, TherapyDB

# query_file = "../tests/data/2025-1-9_civic_failed_pdac_genes.txt"
# <2025-01-13 Mon> The query file below contains genes that had hits in the above
query_file = "/home/shannc/Bio_SDD/chula-stem/tests/data/2025-01-13_pandrugs2_query.txt"
with open(query_file, "r") as f:
    genes = f.read().splitlines()

pandrugs2_cache = (
    "/home/shannc/Bio_SDD/chula-stem/resources/2025-01-13_pandrugs2_cache_pdac.json"
)
P = PanDrugs2(pandrugs2_cache)
P.get_genes(genes, confident=True)

# <2025-01-13 Mon> Now merge the caches
civic_cache = "/home/shannc/Bio_SDD/chula-stem/.cache/civic_cache.json"
civic = pl.read_json(civic_cache).select(TherapyDB.shared_cols)
pd2 = pl.read_json(pandrugs2_cache).select(TherapyDB.shared_cols)
merged = pl.concat([civic, pd2])

all_cache = "/home/shannc/Bio_SDD/chula-stem/.cache/2025-01-13_therapy_data.parquet"
merged.write_parquet(all_cache)

# Schema is
# Schema([('gene', String),
# ('source', String),
# ('disease', List(String)),
# ('therapies', List(String)),
# ('db', String),
# ('db_link', String)])


# * Test get one gene
test = "SMARCA4"
result = P.get_gene(test, confident=True)
