#!/usr/bin/env ipython

import argparse
import shutil
from datetime import date
from pathlib import Path
from typing import Literal, TypeAlias

import tomllib

ON_EXISTS: TypeAlias = Literal["override", "append_date"] | None
TODAY: str = date.today().isoformat()

"""
Helper script for uploading processed data to the CU stem cells NAS

Sample-specific data will be uploaded under the "processed" directory for each
 sample, under a directory for the pipeline name
 EX: For results set "variant_calling_1" on the "PHcase" cohort (PDAC cancer type, exome data)
   The directory containing sample P32's results will be uploaded to the following directory:
       cancer_ngs/exome/PDAC/PHcase/P32/processed/variant_calling_1/

 All other directories/files will be uploaded under the "summary" directory for the specified
 cohort

All details on how to upload are specified by a toml file, see the next
block comment for an example template
"""

"""
samples = ["P1", "P2", "P3"] # (Optional) List of sample directories in results
# These be present in the results directory and will throw an error if they don't
# exist

on_exists = "override" # What to do if a results directory already exists at target path
# Valid options
# - override: override all files
# - append_date: upload instead to <results_name>_<current_date>  
# - null: (default) ignore, do nothing

on_missing = "unassigned" # How to handle samples listed in "samples" that aren't found in the specified cohort
# Valid options:
# - create: create a new directory for the sample name
# - unassigned_create: 
# - unassigned: attempt to look up the sample in the "unassigned" cohort and upload it there
# if found, and do nothing if not found
# - null: (default) ignore, do nothing

[names]
results = "variant_calling_1" # Name of results set
cohort = "PHcase1" # Cohort name
data_modality = "exome"
cancer_type = "PDAC"

[paths]
results = ""    # Path to root directory for results
cancer_ngs = "" # Root path to cancer NGS directory (after mounting)


[sample_mapping] # (Optional) Mapping specifying equivalent samples in the samples
# list and in the cohort i.e. sample aliases
P1 = Patient1
"""


with open("/home/shannc/Bio_SDD/chula-stem/pdac_upload.toml", "rb") as f:
    config = tomllib.load(f)


def validate_config(config: dict) -> None:
    required_keys: dict = {
        "names": ["results", "cohort", "data_modality", "cancer_type"],
        "paths": ["cancer_ngs", "results"],
    }
    for k, v in required_keys.items():
        if not (table := config.get(k)):
            raise ValueError(f"The table {k} is required in the config!")
        for req in v:
            if req not in table:
                raise ValueError(f"The key {req} is required in table {k}!")
            if not table[req]:
                raise ValueError(f"The value to key {req} musn't be null!")


def upload_summaries(
    on_exists: ON_EXISTS, target: Path, source: Path, rname: str, samples: list
) -> None:
    summary_dir = target / "summary" / rname
    if on_exists is not None and on_exists not in {"override", "append_date"}:
        raise ValueError("Invalid value given for `on_exists`")
    if not summary_dir.exists():
        summary_dir.mkdir(parents=True)
    else:
        print(f"Warning: `{rname}` already exists in the cohort summary!")
        if on_exists is None:
            print("No option for `on_exists` specified in config, will not upload")
            return
        elif on_exists == "append_date":
            summary_dir = summary_dir.parent / (f"{rname}_{TODAY}")
    for rdir in source.iterdir():
        if rdir.stem in samples:
            continue
        target_path = summary_dir / rdir.stem
        if not target_path.exists():
            shutil.copytree(rdir, target_path, dirs_exist_ok=False)
        elif on_exists == "override":
            shutil.copytree(rdir, target_path, dirs_exist_ok=True)


def upload_helper(
    on_exists: ON_EXISTS, rname: str, sdir: Path, processed: Path
) -> None:
    """
    Attempt to upload sample directory `sdir` to the `processed` directory of the given sample
    """
    target = processed / rname
    if target.exists() and on_exists == "override":
        shutil.copytree(sdir, target, dirs_exist_ok=True)
    elif target.exists() and on_exists == "append_date":
        target = processed / f"{rname}_{TODAY}"
        shutil.copytree(sdir, target, dirs_exist_ok=False)
    else:
        shutil.copytree(sdir, target, dirs_exist_ok=False)


def upload_samples(
    cohort_dir: Path,
    source: Path,
    rname: str,
    samples: list,
    on_missing: Literal["unassigned", "create", "unassigned_create"] | None,
    on_exists: ON_EXISTS,
) -> None:
    unassigned: Path = cohort_dir.parent.joinpath("unassigned")
    for sample in samples:
        to_upload: Path = source / sample
        if not to_upload.exists():
            raise ValueError(f"Sample {sample} doesn't exist in the results directory")
        target_dir: Path = cohort_dir / sample
        processed: Path = target_dir / "processed"
        if target_dir.exists():
            if not processed.exists():
                processed.mkdir()
        elif on_missing != "create":
            target_dir = unassigned / sample
            processed = target_dir / "processed"
            if not target_dir.exists() and on_missing == "unassigned_create":
                processed.mkdir(parents=True)
            elif not target_dir.exists() and on_missing == "unassigned":
                continue
        else:
            processed.mkdir(parents=True)
        upload_helper(
            on_exists=on_exists, rname=rname, sdir=to_upload, processed=processed
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", required=True)
    config_file = parser.config
    with open(config_file, "rb") as f:
        config: dict = tomllib.load(f)
    validate_config(config)
    paths, names = config["paths"], config["names"]
    rname: str = names["results"]
    source: Path = Path(paths["results"])
    cohort_dir: Path = Path(
        paths["cancer_ngs"]
        / names["data_modality"]
        / names["cancer_type"]
        / names["cohort"]
    )
    samples: list = config.get("samples", [])
    upload_summaries(
        samples=samples, config=config, target=cohort_dir, source=source, rname=rname
    )
