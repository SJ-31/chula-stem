#!/usr/bin/env python3

from pathlib import Path

import skbio.io
from pyhere import here
from skbio import DNA

seq_dir: Path = here("analyses", "data", "tcr_imgt_c")


seqs = [
    next(skbio.io.read(seq_file, "fasta", constructor=DNA, lowercase=True))
    for seq_file in seq_dir.glob("*fasta")
]
for seq in seqs:
    translated = seq.translate()
    header = translated.metadata["id"] + translated.metadata["description"]
    print(f"Header: {header}")
    print(translated)
    print()
    stops = translated.stops()
    if not (sum(stops) == 1 and stops[-1]):
        raise ValueError("Sequence doesn't have terminal stop")
