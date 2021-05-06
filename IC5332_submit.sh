#!/bin/bash
#
#SBATCH --job-name=IC5332_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -p OBS,PREEMPTION
#SBATCH --mem-per-cpu=1000
#SBATCH --output=IC5332_output.out
#SBATCH --error=IC5332_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

srun python fitmodel.py IC5332
