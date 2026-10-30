#!/usr/bin/env ipython

import json
from collections import defaultdict
from pathlib import Path
from typing import Literal

import duckdb
import polars as pl
from attrs import Factory, define, field, validators

# [2026-10-30 Fri]
# Plan for using maveDB data (and any other MAVE data)
# Just add the data to all possible variants in the report, showing
# the relevant experiments and the MAVE score
# If the variant is of uncertain significance, make a note of it


@define
class MaveDB:
    dir: Path = field(converter=Path)
    mapped: Path = field(
        init=False, default=Factory(lambda x: x / "mapped", takes_self=True)
    )
    csv: Path = field(init=False, default=Factory(lambda x: x / "csv", takes_self=True))
    va: Path = field(init=False, default=Factory(lambda x: x / "va", takes_self=True))
    main: dict = field(
        init=False,
        default=Factory(
            lambda x: json.load(open(x / "main.json", "rb")), takes_self=True
        ),
    )
    score_sets: pl.DataFrame = field()

    @score_sets.default
    def _make_score_sets(self) -> pl.DataFrame:
        # TODO: flatten the main dictionary into a df of metadata for
        # score sets, to facilitate easy lookups
        tmp = {
            "urn": [],
            "targetGenes": [],
            "targetGenes.HGNC": [],
            "targetGenes.UniprotID": [],
        }
        external_ids = defaultdict(dict)
        for exp_set in self.main["experimentSets"]:
            for exp in exp_set["experiments"]:
                for score_set in exp["scoreSets"]:
                    for tg in score_set["targetGenes"]:
                        name = tg["name"]
                        tmp["targetGenes"].append(name)
                        tmp["targetGenes.UniprotID"].append(
                            tg["uniprotIdFromMappedMetadata"]
                        )
                        tmp["targetGenes.HGNC"].append(tg["mappedHgncName"])
                        tmp["urn"].append(score_set["urn"])
                        for ext in tg["externalIdentifiers"]:
                            db = ext["identifier"]["dbName"]
                            id = ext["identifier"]["identifier"]
                            external_ids[db][id] = name
        df = pl.DataFrame(tmp)
        for k, v in external_ids.items():
            ext_df = pl.DataFrame(
                {"targetGenes": v.values(), f"targetGenes.{k}": v.keys()}
            )
            df = df.join(ext_df, on="targetGenes", how="left")

        return df

    @dir.validator
    def _check_dump(self, _attr, val: Path):
        assert val.exists()
        assert val.is_dir()
        assert (val / "main.json").exists()
        assert (val / "va").exists()
        assert (val / "csv").exists()

    def __attrs_post_init__(self):
        pass

    def __repr__(self) -> str:
        return repr(self.dir)

    def __truediv__(self, path) -> Path:
        return self.dir / path

    def urn2file(self, urn: str, which: Literal["va", "mapped", "csv"]) -> Path:
        name = urn.replace(":", "-")
        if which == "mapped":
            suffix = ".mapped-variants.json"
        elif which == "mapped":
            suffix = ".va.ndjson"
        else:
            suffix = ".scores.csv"
        lookup = self / which / name / suffix
        if lookup.exists():
            return lookup
        raise FileNotFoundError(f"{which} file for URN {urn} is missing")

    # TODO: read the VRS documentation

    @staticmethod
    def parse_mapping(map_file: Path) -> pl.DataFrame:
        tmp = {}
        with open(map_file, "rb") as f:
            txt = json.load(f)

    def gene_score_set(self, gene: str, col: str = "targetGenes.HGNC") -> pl.DataFrame:
        urns = self.score_sets.filter(pl.col(col) == gene)["urn"]
        for urn in urns:
            scores = pl.read_csv(self.urn2file(urn, "csv"))
            try:
                va_text = self.urn2file(urn, "")
                with open(self.urn2file(urn, "va")) as va_file:
                    va_tmp = json.load(va_file)
            except FileNotFoundError:
                pass

        pass


M = MaveDB(
    "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/reference/variants/mavedb-dump.2026062418131"
)
