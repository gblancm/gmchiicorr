#!/bin/bash
#
#SBATCH --job-name=NGC3627_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=168:00:00
#SBATCH -p any
#SBATCH --mem-per-cpu=2000
#SBATCH --output=NGC3627_output.out
#SBATCH --error=NGC3627_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load python3/3.8.5
srun python3 fitmodel.py NGC3627
