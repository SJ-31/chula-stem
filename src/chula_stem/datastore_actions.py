import datetime
import os
import re
import shutil
from collections.abc import Sequence
from functools import reduce
from pathlib import Path
from typing import Any

import click
import numpy as np
import pandas as pd
import yaml


def get_confirmation(
    prompt, valid_responses=("y", "n"), mapping: None | dict = None
) -> Any:
    if not mapping:
        mapping = {"y": True, "n": False}
    while (response := input(f"{prompt}: ").lower()) not in valid_responses:
        print(f"Invalid response. Available options: {valid_responses}")
    return mapping[response]


def dir_is_empty(dir: Path) -> bool:
    if not dir.exists():
        return True
    return next(dir.iterdir(), None) is None


def find_cases(tumor_dir: Path, modality: str, date: str) -> pd.DataFrame:
    data = {
        "cohort": [],
        "case_name": [],
        "path": [],
        "tumor_type": [],
        "modality": [],
        "date_received": [],
        "has_pbmc": [],
        "has_tumor": [],
        "has_raw": [],
        "has_processed": [],
    }
    for cohort in tumor_dir.iterdir():
        for case_dir in cohort.iterdir():
            if case_dir.stem == "summary":
                continue
            for name in ("processed", "raw"):
                dir = case_dir.joinpath(name)
                is_empty = dir_is_empty(dir)
                if dir.exists() and not is_empty:
                    data[f"has_{name}"].append("T")
                else:
                    data[f"has_{name}"].append("F")
                if name == "raw" and not is_empty:
                    contents = list(dir.iterdir())
                    for suffix, stype in zip(("_B", "_C"), ("pbmc", "tumor")):
                        check = list(map(lambda x: suffix in x.stem, contents))
                        if any(check):
                            data[f"has_{stype}"].append("T")
                        else:
                            data[f"has_{stype}"].append("F")
                elif name == "raw":
                    data["has_pbmc"].append("F")
                    data["has_tumor"].append("F")
            data["modality"].append(modality)
            data["case_name"].append(case_dir.stem)
            data["cohort"].append(cohort.stem)
            data["tumor_type"].append(tumor_dir.stem)
            data["path"].append(str(case_dir.absolute()))
            data["date_received"].append(date)
    return pd.DataFrame(data)


def get_dirs(root: Path) -> Sequence[Path]:
    lst = []
    for m in root.iterdir():
        if not os.access(m, os.R_OK) or not m.is_dir() or m.stem == "TMP":
            continue
        lst.append(m)
    return lst


def get_entry_df(root: Path, date: str) -> pd.DataFrame:
    dfs: list = []
    for modality in get_dirs(root):
        for tumor_type in modality.iterdir():
            dfs.append(
                find_cases(tumor_dir=tumor_type, modality=modality.stem, date=date)
            )
    return pd.concat(dfs)


def parse_date(date_string) -> datetime.datetime | None:
    """Attempt to parse `date_string` into a datetime object (ISO8601),
    returning None upon failure"""
    try:
        return datetime.datetime.strptime(date_string, "%Y-%m-%d")
    except ValueError:
        return None


def date_x_newer(x: str, y: str) -> bool:
    """True if date string x is more recent than date string y, False
    otherwise
    """
    parsed_x, parsed_y = parse_date(x), parse_date(y)
    if parsed_x is not None and parsed_y is None:
        return True
    elif parsed_x is not None and parsed_y is not None and (parsed_x > parsed_y):
        return True
    return False


def write_entry_df(root, date, previous_path: str | None, output: str) -> None:
    df: pd.DataFrame = get_entry_df(root=root, date=date)
    if previous_path is not None:
        previous: pd.DataFrame = pd.read_csv(previous_path)
        previous.loc[:, "version"] = "old"

        for frame in [previous, df]:
            id_cols = [
                frame[x] for x in ["cohort", "case_name", "tumor_type", "modality"]
            ]
            primary_key = reduce(lambda x, y: x.str.cat(y), id_cols)
            frame.loc[:, "key"] = primary_key

        new_entries: pd.DataFrame = df.loc[~df["key"].isin(previous["key"]), :].drop(
            columns="key"
        )
        # Permit newer versions if...
        # - The previous entry had a empty date, but the new version doesn't
        # - Both entries have dates and the new version is more recent
        joined = pd.merge(
            previous, df, how="inner", on="key", suffixes=("_prev", "_new")
        )
        df = df.loc[df["key"].isin(previous["key"]), :]

        get_from_new = np.array(
            map(
                lambda x: date_x_newer(x[0], x[1]),
                zip(joined["date_received_new"], joined["date_received_prev"]),
            )
        )
        newer_dates = df.loc[get_from_new, :].drop(columns="key")
        kept_dates = previous.loc[~get_from_new, :].drop(columns="key")

        df = pd.concat([new_entries, newer_dates, kept_dates])
    df.to_csv(output, index=False)


def store_away(config_path: str, sink_root, interactive, noconfirm, output):
    """Store directories into NGS data root by name

    Parameters
    ----------
    config : YAML file providing parameters and options for storing away a set of files
        e.g. pipeline results. Config files should be written to be specific for a
        given cohort, modality and tumor type


    Notes
    -----
    config support the following keys:
    - source_root : str
        Required key which specifies directory containing data to store away.
        Subdirectories must be case names
    - copy : bool
        if True, copy files from ``source_root`` instead of renaming them into
        the sink
    - modality : str
        Required key to specify modality directory in ``sink_root``,
    - condition : str
        Specify condition of cases e.g. a tumor type. Default = "unknown"
    - cohort : str
        Specify case cohort. Default = "unassigned"
    - case_mappings : map
        Maps specific case names to new names for their destination in ``sink_root``
    - case_store_dir : str
        Subdirectory in case to put data in. Either "processed" or "raw"
    - ignored : boolean
        list of case subdirectories to ignore
    - case_regexp : str
        Regexp used to check if a subdirectory is a case and should be moved
    - store_name : str
        Subdirectory under ``case_store_dir`` to move data into. If unspecified, contents
        are placed with no directory
    - store_name_is_prefix : bool
        If True, the subdirectory for storing data is then "{store_name}-{case_name}"
    """
    with open(config_path, "r") as y:
        config: dict = yaml.safe_load(y)

    sink_root_path = Path(sink_root)
    store_path: Path = Path(sink_root)
    command_log: dict = {"old": [], "new": []}

    # Create or validate datastore directories to copy/move into
    for subdir_name, default in zip(
        ["modality", "condition", "cohort"], [None, "unknown", "unassigned"]
    ):
        given = config.get(subdir_name, default)
        if not given:
            raise ValueError("Modality must be specified in config!")
        available = get_dirs(store_path)
        store_path = store_path.joinpath(given)
        if given not in available and not Path(store_path).exists():
            print(f"{subdir_name} {given} not found in {store_path.parent}")
            if not interactive and not noconfirm:
                return
            if noconfirm or get_confirmation(f"Create directory {store_path}? Y/N"):
                store_path.mkdir()

    source_root = Path(config["source_root"])
    if not source_root.exists():
        raise ValueError("source_root not found")

    ignored = config.get("ignored", set())
    do_copy = config.get("copy", False)
    case_regexp = config.get("case_regexp", "")
    store_subdir = config.get("case_store_dir", "processed")
    case_mappings: dict = config.get("case_mappings", {})

    # Copy/move cases into datastore
    for case in get_dirs(source_root):
        if case.stem in ignored:
            continue
        elif case_regexp and not re.match(case_regexp, case.stem):
            print(f"Subdir {case} did not match case regexp, ignoring...\n")
            continue
        case_alias = case_mappings.get(case.stem, case.stem)
        cur_path = store_path.joinpath(case_alias)
        if not cur_path.exists():
            if not interactive and not noconfirm:
                raise ValueError(
                    "Cannot create new case directories if not interactive or create = False"
                )
            print(f"WARNING: directory for case {case_alias} does not exist")
            if noconfirm or get_confirmation(f"Create {cur_path}? Y/N"):
                cur_path.mkdir()
                print()
            else:
                continue
        cur_path = cur_path.joinpath(store_subdir)  # e.g. case/raw
        cur_path.mkdir(exist_ok=True)
        if store_name := config.get("store_name"):
            if config.get("store_name_is_prefix", False):
                store_name = f"{store_name}-{case_alias}"
            cur_path = cur_path.joinpath(store_name)
            cur_path.mkdir(exist_ok=True)
        for file in case.iterdir():
            rel = file.relative_to(source_root)
            exists = cur_path.joinpath(file.name).exists()
            if exists and not interactive and not noconfirm:
                print(f"WARNING: file {rel} already exists in {cur_path}. Skipping...")
                continue
            elif exists and interactive:
                if not get_confirmation(
                    f"WARNING: file {rel} already exists in {cur_path}. [s]kip or add with [u]nique name?",
                    valid_responses=("s", "u"),
                    mapping={"u": True, "s": False},
                ):
                    print()
                    continue
            if not do_copy and not exists:
                shutil.move(file, cur_path)
            elif not do_copy and exists:
                shutil.move(file, make_unique_file(file.name, cur_path))
            elif file.is_dir():
                shutil.copytree(file, cur_path, dirs_exist_ok=True)
            else:
                shutil.copy(file, cur_path)
        command_log["old"].append(f"{source_root}/{case}")
        command_log["new"].append(cur_path.relative_to(sink_root_path))
    df = pd.DataFrame(command_log)
    df.to_csv(output, index=False)


def make_unique_file(file: str, dir: Path) -> Path:
    current = dir.joinpath(file)
    i: int = 1
    while current.exists():
        suffixes = current.suffixes
        if len(suffixes) > 0:
            base = reduce(lambda x, y: x.replace(y, ""), suffixes, current.name)
            new_name = "".join([base, f"_{i}"] + suffixes)
        else:
            new_name = f"{current.name}_{i}"
        current = dir.joinpath(new_name)
        i += 1
    return current


@click.command()
@click.option(
    "-d", "--date", required=False, help="Entry for the 'date' column", default=""
)
@click.option("-r", "--root", required=False, help="Root of NGS data", default=".")
@click.option(
    "-s",
    "--store_config",
    required=False,
    help="YAML file specifying file storage behavior. Required for the `store_away` action",
    default=None,
)
@click.option(
    "-a",
    "--action",
    required=True,
    help="Action to perform",
    default="record_cases",
    type=click.Choice(["record_cases", "store_away"]),
)
@click.option(
    "-i",
    "--interactive",
    required=False,
    help="Whether to operate interactively",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-n",
    "--noconfirm",
    required=False,
    help="""
    Do not require user confirmation for certain actions.
    Specifically, non-existent directories specified in the store_config will
    be created automatically, and existing files in the datastore will be overwritten
    if new versions are present in the source directory
    """,
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-o", "--output", required=False, help="File to write record to", default=None
)
@click.option(
    "-p",
    "--previous_record",
    required=False,
    help="Previously generated sample manifest. If supplied, the new manifest will only override data in old entries if the date string has updated",
    default=None,
)
def main(
    date, root, action, previous_record, output, store_config, interactive, noconfirm
) -> None:
    if not output:
        d = datetime.date.today().isoformat()
        output = f"{d}-{action}-autogen.csv"
    if action == "record_cases":
        write_entry_df(
            root=Path(root), date=date, previous_path=previous_record, output=output
        )
    elif action == "store_away":
        if store_config is None:
            raise ValueError("--store_config (-s) must be provided for this action!")
        store_away(
            config_path=store_config,
            sink_root=root,
            interactive=interactive,
            noconfirm=noconfirm,
            output=output,
        )


if __name__ == "__main__":
    main()
