#!/bin/bash
#SBATCH --job-name=vertebra
#SBATCH --account=pi-bradnelson
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err
#SBATCH --time=01:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=200G

module load python

conda init --all
source /home/yuanluo/.bashrc
conda activate /home/yuanluo/.conda/envs/bats
python vertebra32.py