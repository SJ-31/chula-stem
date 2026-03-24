#!/usr/bin/env ipython

from pathlib import Path

import yte
from pyhere import here

env = yte.process_yaml(here("analyses", "cluster_reproducibility", "env.yaml"))
workdir = here(*env["workdir"])
outdir = here(*env["outdir"])
