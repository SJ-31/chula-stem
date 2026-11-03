#!/usr/bin/env ipython

import json
from collections import defaultdict
from pathlib import Path
from typing import Literal

import polars as pl
from attrs import Factory, define, field

# [2026-10-30 Fri]
# Plan for using maveDB data (and any other MAVE data)
# Just add the data to all possible variants in the report, showing
# the relevant experiments and the MAVE score
# If the variant is of uncertain significance, make a note of it

# https://www.mavedb.org/docs/mavedb/finding-data/downloading.html#archive-contents


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
        tmp = {
            "urn": [],
            "title": [],
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
                        tmp["title"].append(score_set["title"])
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
        elif which == "va":
            suffix = ".va.ndjson"
        else:
            suffix = ".scores.csv"
        return self / which / f"{name}{suffix}"

    @staticmethod
    def get_va(va_file: Path) -> pl.DataFrame:
        """ """
        tmp = {
            "variant_urn": [],
            "va_spec": [],
            "variant_ga4gh": [],
            "evidence_direction": [],
            "classification": [],
        }
        for json_line in va_file.read_text().splitlines():
            data = json.loads(json_line)
            anno = data["annotation"]
            if not anno:
                continue
            match anno:
                case {"proposition": prop}:
                    match prop:
                        case {"subjectVariant": subject} | {"subject": subject}:
                            tmp["variant_ga4gh"].append(subject["id"])
                    tmp["va_spec"].append(anno["proposition"]["type"])
                    tmp["evidence_direction"].append(anno["direction"])
                    tmp["classification"].append(
                        anno["classification"]["primaryCoding"]["code"]
                    )
                case {"focusVariant": focus} | {"focus": focus}:
                    if "id" not in focus:
                        continue
                    tmp["variant_ga4gh"].append(focus["id"])
                    tmp["va_spec"].append(anno["type"])
                    tmp["classification"].append(None)
                    tmp["evidence_direction"].append(None)
            tmp["variant_urn"].append(data["variant_urn"])
        return pl.DataFrame(tmp)

    @staticmethod
    def get_mapping(map_file: Path) -> pl.DataFrame:
        """
        Retrieve successfully mapped variants from `map_file`
        """
        rows = {"variant_urn": [], "current": []}
        with open(map_file, "rb") as f:
            data: dict = json.load(f)

        tmp_exprs = defaultdict(list)

        for rec in data:
            post = rec.get("postMapped")
            pre = rec.get("preMapped")
            if rec["errorMessage"] or (not post and not pre):
                continue
            exprs = post.get("expressions") or []
            vu = rec.get("variantUrn")
            if not vu:
                continue
            rows["variant_urn"].append(vu)
            rows["current"].append(rec.get("current"))
            for expr in exprs:
                tmp_exprs["variant_urn"].append(vu)
                tmp_exprs[expr["syntax"]].append(expr["value"])

        expr_df = pl.DataFrame(tmp_exprs)
        df = pl.DataFrame(rows).join(expr_df, on="variant_urn")
        return df

    def gene_score_set(self, gene: str, col: str = "targetGenes.HGNC") -> pl.DataFrame:
        rel = self.score_sets.filter(pl.col(col) == gene)
        dfs = []
        for urn in rel["urn"]:
            scores = pl.read_csv(self.urn2file(urn, "csv"))
            if (mapped_file := self.urn2file(urn, "mapped")).exists():
                mapped = self.get_mapping(mapped_file)
                scores = scores.join(
                    mapped, left_on="accession", right_on="variant_urn"
                )
            if (va_file := self.urn2file(urn, "va")).exists():
                va = self.get_va(va_file)
                scores = scores.join(va, left_on="accession", right_on="variant_urn")
            scores = scores.with_columns(pl.lit(urn).alias("urn"))
            dfs.append(scores)
        final = pl.concat(dfs, how="diagonal_relaxed").join(rel, on="urn")
        return final


# TODO: when formatting the report, create a link to the variant on the
# mavedb website from the ga4gh digest
# And the source of the variant can also be the main score set URN
# e.g. https://mavedb.org/score-sets/urn:mavedb:00000001-a-1
