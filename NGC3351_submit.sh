#!/bin/bash
#
#SBATCH --job-name=gmc_hii_corr_mod
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24 
#SBATCH --time=24:00:00     
#SBATCH -p OBS,SHARED
#SBATCH --mem-per-cpu=4000
#SBATCH --output=NGC3351_output.out 
#SBATCH --error=NGC3351_error.out 
#SBATCH --mail-user=gblancm@carnegiescience.edu 
#SBATCH --mail-type=END,FAIL 

#srun module load python
srun python gmc_hii_corr_model.py NGC3351