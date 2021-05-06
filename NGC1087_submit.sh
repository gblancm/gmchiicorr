#!/bin/bash
#
#SBATCH --job-name=NGC1087_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -p OBS
#SBATCH --mem-per-cpu=4000
#SBATCH --output=NGC1087_output.out
#SBATCH --error=NGC1087_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

srun python fitmodel.py NGC1087
