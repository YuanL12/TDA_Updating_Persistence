#!/bin/sh

#SBATCH --job-name=levelset-3d
#SBATCH --account=pi-bradnelson
#SBATCH --output=levelset_3d.out
#SBATCH --error=levelset_3d.err
#SBATCH --time=10:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G

module load python
conda init --all
conda activate bats

# use time to get memory use information
# /usr/bin/time -v python levelset_3d.py
python levelset_3d.py
