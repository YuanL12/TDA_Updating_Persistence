#!/bin/bash
#SBATCH --job-name=vertebra64
#SBATCH --account=pi-bradnelson
#SBATCH --output=vertebra64.out
#SBATCH --error=vertebra64.err
#SBATCH --time=1:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G

module load python

conda init --all
conda activate tda

python ../vertebra64.py
