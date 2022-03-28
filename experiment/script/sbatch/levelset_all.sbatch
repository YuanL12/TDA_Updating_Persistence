#!/bin/sh

#SBATCH --job-name=levelset-all
#SBATCH --account=pi-bradnelson
#SBATCH --output=levelset_all.out
#SBATCH --error=levelset_all.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G

module load python
conda init --all
conda activate tda

# use time to get memory use information
# /usr/bin/time -v python levelset_3d.py
python ../levelset_all.py
