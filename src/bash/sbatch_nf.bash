#!/usr/bin/env bash

d=$(date +'%A_%b-%d-%Y')
sbatch -e "sbatch_nf_stderr_${d}.txt" -o "sbatch_nf_stdout_${d}.txt" \
    "$HOME/chula-stem/src/bash/nf_slurm.bash" "$@"
