#!/usr/bin/env python

import argparse
import csv
import pprint
import shutil
from collections import Counter, defaultdict
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
summary_ignore = [""] # (Optional) directory/file names to ignore in summary

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
cancer_type = "pdac" # Cancer type acronym (should be lowercase)

[paths]
results = ""    # Path to root directory for results
cancer_ngs = "" # Root path to cancer NGS directory (after mounting)

[sample_mapping] # (Optional) Mapping specifying equivalent samples in the samples
# list and in the cohort i.e. sample aliases
P1 = Patient1

[cohort_override] # (Optional) Mapping specifying alternate cohort for a specific sample
# P1 = "unassigned"
# CEN2-P3 = "DCEN2"
"""


def validate_config(config: dict, dry_run: bool) -> tuple[Path, Path, list, dict]:
    """Validate config file for correct entries

    Returns
    -------
    Tuple of (source_path, cohort_dir, sample_list)
    """
    required_keys: dict = {
        "names": ["results", "cohort", "data_modality", "cancer_type"],
        "paths": ["cancer_ngs", "results"],
    }
    for k, v in required_keys.items():
        if not (table := config.get(k)):
            raise ValueError(f"The table {k} is required in the config")
        for req in v:
            if req not in table:
                raise ValueError(f"The key {req} is required in table {k}")
            if not table[req]:
                raise ValueError(f"The value to key {req} musn't be null")
    paths, names = config["paths"], config["names"]
    source = Path(paths["results"])
    if not source.exists():
        raise ValueError(f"The results source {source} doesn't exist")
    cohort_dir: Path = Path(paths["cancer_ngs"])
    for key in ["data_modality", "cancer_type", "cohort"]:
        if not (try_path := cohort_dir.joinpath(names[key])).exists():
            print(f"The path {try_path} given by key {key} doesn't exist")
            print(f"\tCreating path {try_path}")
            smart_mkdir(try_path, dry_run=dry_run)
        cohort_dir = try_path

    sample_mapping = config.get("sample_mapping", {})
    cohort_override: dict = config.get("cohort_override", {})
    samples: list = config.get("samples", [])
    default_cohort_name = names["cohort"]
    samples_check: dict = defaultdict(list)
    for s in samples:
        cur_cohort = cohort_override.get(s, default_cohort_name)
        if not source.joinpath(s).exists():
            raise ValueError(f"Sample {s} doesn't exist in the results directory")
        if s in sample_mapping:
            print(f"Using remapping {s} = {sample_mapping[s]}")
            s = sample_mapping[s]
        samples_check[cur_cohort].append(s)
    for k, v in samples_check.items():
        if len(v) != len(set(v)):
            message = f"Sample mapping in cohort {k} causes name conflicts"
            message += (
                f"\n\tDuplicated names: {[k for k,v in Counter(v).items() if v> 1]}"
            )
            raise ValueError(message)
    print("\nFinal cohort assignment:")
    pprint.pprint(samples_check)
    print()
    return source, cohort_dir, samples, sample_mapping


def smart_mkdir(dir: Path, parents=False, dry_run: bool = False) -> bool:
    if not dry_run and not dir.exists():
        dir.mkdir(parents=parents)
        return True
    return False


def upload_summaries(
    on_exists: ON_EXISTS,
    target: Path,
    source: Path,
    rname: str,
    config: dict,
    dry_run: bool,
) -> list[dict]:
    samples = config.get("samples", [])
    to_ignore = set(samples) | set(config.get("summary_ignore", []))
    summary_dir = target / "summary" / rname
    tracker: list[dict] = []
    if on_exists is not None and on_exists not in {"override", "append_date"}:
        raise ValueError("Invalid value given for `on_exists`")
    if not smart_mkdir(summary_dir, parents=True):
        print(f"Warning: `{rname}` already exists in the cohort summary!")
        if on_exists is None:
            print("\tNo option for `on_exists` specified in config, will not upload")
            return []
        elif on_exists == "append_date":
            summary_dir = summary_dir.parent / (f"{rname}_{TODAY}")
    for rdir in source.iterdir():
        if rdir.stem in to_ignore:
            continue
        target_path = summary_dir / rdir.stem
        template = {"old": rdir, "status": "success", "new": target_path}
        if target_path.exists() and on_exists == "override":
            print(f"WARNING: overriding {target_path}")
            template["override"] = True
        else:
            template["override"] = False
        if target_path.exists() and not on_exists:
            template["status"] = "failure"
            template["new"] = None
        elif dry_run:
            print(f"{rdir}->{target_path}")
        elif rdir.is_file():
            shutil.copy2(rdir, target_path)
        else:
            shutil.copytree(rdir, target_path, dirs_exist_ok=True)
        tracker.append(template)
    return tracker


def upload_helper(
    on_exists: ON_EXISTS, rname: str, sdir: Path, processed: Path, dry_run: bool
) -> tuple[Path, bool]:
    """
    Attempt to upload sample directory `sdir` to the `processed` directory of the given sample

    Returns
    -------
    Tuple of (destination, boolean which is True if destination files were overriden)
    """
    target = processed / rname
    if target.exists() and on_exists == "append_date":
        target = processed / f"{rname}_{TODAY}"
        if dry_run:
            print(f"{sdir}->{target}")
        else:
            shutil.copytree(sdir, target, dirs_exist_ok=False)
    elif dry_run:
        print(f"{sdir}->{target}")
    else:
        shutil.copytree(sdir, target, dirs_exist_ok=True)
    return target, on_exists == "override" and target.exists()


def upload_samples(
    initial_cohort_dir: Path,
    source: Path,
    rname: str,
    samples: list,
    sample_mapping: dict,
    cohort_override: dict,
    on_missing: Literal["unassigned", "create", "unassigned_create"] | None,
    on_exists: ON_EXISTS,
    dry_run: bool,
) -> list:
    unassigned: Path = initial_cohort_dir.parent.joinpath("unassigned")
    sample_tracker: list[dict] = []
    print("---")
    for sample in samples:
        to_upload: Path = source / sample
        template = {"sample": sample, "old": to_upload}
        should_upload: bool = True
        if c_override := cohort_override.get(sample):
            print("Overriding cohort...")
            cohort_dir = initial_cohort_dir.parent.joinpath(c_override)
        else:
            cohort_dir = initial_cohort_dir
        sample = sample_mapping.get(sample, sample)
        target_dir: Path = cohort_dir / sample
        processed: Path = target_dir / "processed"
        if not target_dir.exists():
            print(f"WARNING: sample {sample} not found in {cohort_dir}")

        if target_dir.exists():
            smart_mkdir(processed, dry_run=dry_run)
        elif not on_missing:
            print("\t`on_missing` behavior not specified, skipping")
            should_upload = False
        elif on_missing != "create":
            target_dir = unassigned / sample
            processed = target_dir / "processed"
            if not target_dir.exists() and on_missing == "unassigned_create":
                print("\tNot found in `unassigned`, creating sample directory")
                smart_mkdir(processed, parents=True, dry_run=dry_run)
            elif not target_dir.exists() and on_missing == "unassigned":
                print("\tNot found in `unassigned`, skipping")
                should_upload = False
        elif on_missing == "create":
            print("\tCreating sample directory")
            smart_mkdir(processed, parents=True, dry_run=dry_run)
        uploaded, override = (
            upload_helper(
                on_exists=on_exists,
                rname=rname,
                sdir=to_upload,
                processed=processed,
                dry_run=dry_run,
            )
            if should_upload
            else (None, False)
        )
        template["status"] = "failure" if not uploaded else "success"
        template["override"] = override
        template["new"] = uploaded
        sample_tracker.append(template)
        print("---")
    return sample_tracker


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument(
        "-l", "--log", required=False, default=f"upload_record-{TODAY}.csv"
    )
    parser.add_argument(
        "-d",
        "--dry_run",
        required=False,
        default=False,
        action="store_true",
    )
    args = parser.parse_args()
    config_file = args.config

    with open(config_file, "rb") as f:
        config: dict = tomllib.load(f)
    source: Path
    cohort_dir: Path
    source, cohort_dir, samples, sample_mapping = validate_config(
        config, dry_run=args.dry_run
    )
    print(f"Default cohort: {cohort_dir}")

    rname: str = config["names"]["results"]

    on_exists: ON_EXISTS = config.get("on_exists")
    summaries_tracker = [
        dict({"sample": None}, **r)
        for r in upload_summaries(
            config=config,
            target=cohort_dir,
            source=source,
            rname=rname,
            on_exists=on_exists,
            dry_run=args.dry_run,
        )
    ]
    samples_tracker = upload_samples(
        samples=samples,
        source=source,
        rname=rname,
        initial_cohort_dir=cohort_dir,
        cohort_override=config.get("cohort_override", {}),
        on_exists=on_exists,
        on_missing=config.get("on_missing"),
        sample_mapping=sample_mapping,
        dry_run=args.dry_run,
    )

    with open(args.log, "w") as logfile:
        fields = ["sample", "old", "new", "status", "override"]
        writer = csv.DictWriter(logfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries_tracker)
        writer.writerows(samples_tracker)
