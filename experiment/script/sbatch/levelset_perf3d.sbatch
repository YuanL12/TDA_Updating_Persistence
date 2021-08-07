#!/bin/sh

#SBATCH --job-name=levelset-perf3d
#SBATCH --account=pi-bradnelson
#SBATCH --output=levelset_perf3d.out
#SBATCH --error=levelset_perf3d.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64G

module load python
conda init --all
conda activate bats

# use time to get memory use information
# /usr/bin/time -v python levelset_3d.py
python levelset_perf3d.py
