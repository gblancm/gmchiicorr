#!/bin/bash
#
#SBATCH --job-name=NGC1385_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=168:00:00
#SBATCH -p any
#SBATCH --mem-per-cpu=2000
#SBATCH --output=NGC1385_output.out
#SBATCH --error=NGC1385_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

srun python3 fitmodel.py NGC1385
