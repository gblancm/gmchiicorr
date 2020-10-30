#!/bin/bash
#
#SBATCH --job-name=gmc_hii_corr_mod
#SBATCH --output=output.txt
#
#SBATCH --ntasks=5
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100

srun module load python
srun python gmc_hii_corr_model.py 
