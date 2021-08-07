#!/bin/sh

#SBATCH --job-name=levelset-3d-cubical
#SBATCH --account=pi-bradnelson
#SBATCH --output=levelset_3d_cubical.out
#SBATCH --error=levelset_3d_cubical.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64G

module load python
conda init --all
conda activate bats

# use time to get memory use information
# /usr/bin/time -v python levelset_3d_cubical.py
python levelset_3d_cubical.py
