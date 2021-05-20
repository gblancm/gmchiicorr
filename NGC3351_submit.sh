#!/bin/bash
#
#SBATCH --job-name=NGC3351_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -p OBS,PREEMPTION,SHARED
#SBATCH --mem-per-cpu=2000
#SBATCH --output=NGC3351_output.out
#SBATCH --error=NGC3351_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

srun python fitmodel.py NGC3351
