#!/bin/bash
#
#SBATCH --job-name=NGC3627_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -p OBS,SHARED
#SBATCH --mem-per-cpu=4000
#SBATCH --output=NGC3627_output.out
#SBATCH --error=NGC3627_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=END,FAIL

srun python fitmodel.py NGC3627
