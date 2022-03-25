#!/bin/sh

#SBATCH --job-name=levelset-mnist
#SBATCH --account=pi-bradnelson
#SBATCH --output=levelset_mnist.out
#SBATCH --error=levelset_mnist.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G

module load python
conda init --all
conda activate tda

# use time to get memory use information
# /usr/bin/time -v python levelset_3d.py
python ../levelset_mnist.py
