#!/usr/bin/env ipython

import polars as pl

# <2025-01-07 Tue> Want to try to improve caching
civic_json = "/home/shannc/Bio_SDD/chula-stem/tests/civic.json"
civic = pl.read_json(civic_json)
