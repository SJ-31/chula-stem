#!/usr/bin/env ipython

from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import mudata as md
import numpy as np
import polars as pl
import polars.selectors as cs
import pydna.tm as tm
from pydna.design import Amplicon, pcr, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer


def calc_tm(
    seq: str,
    how: Literal["neb", "default"] = "default",
    prodcode: str = "hstaq-0",
    **kws,
):
    """Wrapper function around pydna tm calculations
    See https://tmapi.neb.com/docs/productcodes for product codes
    """
    if how == "default":
        return tm.tm_default(seq, **kws)
    elif how == "neb":
        return tm.tm_neb(seq, prodcode=prodcode, **kws)


def get_cdr3_fusion_primers(
    airr_data: dict,
    chain: str,
    how: Literal["pydna", "manual"] = "pydna",
    forward_len: int = 30,
    reverse_len: int = 30,
    tm_func: Callable[[str], float] = tm.tm_default,
    target_tm: float = 55.0,
) -> dict:
    full: Dseqrecord = Dseqrecord(airr_data[f"{chain}_sequence"])
    p1_start, p1_end = (
        airr_data[f"{chain}_v_sequence_start"],
        airr_data[f"{chain}_cdr3_end"],
    )
    p2_start, p2_end = (
        airr_data[f"{chain}_cdr3_start"],
        airr_data[f"{chain}_j_sequence_end"],
    )
    results = {}
    for i, (start, end) in enumerate(zip((p1_start, p1_end), (p2_start, p2_end))):
        if how == "pydna":
            target = full[start:end]
            amp = primer_design(target, target_tm=target_tm, tm_func=tm_func)
        elif how == "manual":
            fwd = full[start : start + forward_len]
            rev = full[end - reverse_len : end].reverse_complement()
            primers = (Primer(fwd, id=f"f{end-start}"), Primer(rev, id=f"r{end-start}"))
            amp = pcr(primers, full)
        results[f"pair {i+1}"] = amp

    return results


def write_primer_report(
    primers: dict[str, Amplicon], savepath: Path | None = None, **kws
) -> None:
    lines = []
    for name, amp in primers.items():
        lines.extend([f"Primer set: {name}"])
        lines.extend([amp.figure(), ""])
        for p, name in zip(amp.primers(), ("Forward", "Reverse")):
            if p is not None:
                lines.append(name)
                seq = str(p.seq)
                tm = calc_tm(seq, **kws)
                lines.append(f"\tSequence: {seq}")
                lines.append(f"\tLength: {len(seq)}\tTM: {tm}\tGC content: {p.gc()}")
        lines.append("\n")
    text = "\n".join(lines)
    if savepath is None:
        print(text)
    else:
        savepath.write_text(text)
