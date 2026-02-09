#!/bin/bash
#SBATCH --job-name=marimo
#SBATCH --qos=cpu24h
#SBATCH --mem=20G

echo "Connect with: ssh -L 8023:localhost:8023 $USER@$(hostname -f) ssh -L 8023:localhost:8023 $(hostname)"
marimo edit notebook.py --port 8023 --headless --no-token

# ssh -L 8045:localhost:8023 shannc@161.200.107.77 \
#     ssh -L 8023:localhost:8023 cuaim-fat-01.cuaim.hpc
