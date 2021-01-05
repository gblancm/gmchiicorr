#!/bin/bash
#
#SBATCH --job-name=gmc_hii_corr_mod
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=24  
#SBATCH --time=24:00:00     
#SBATCH -p OBS,SHARED
#SBATCH --mem-per-cpu=2000
#SBATCH --output=output.out 
#SBATCH --error=error.out 
#SBATCH --mail-user=gblancm@carnegiescience.edu 
#SBATCH --mail-type=END,FAIL 

#srun module load python
srun python gmc_hii_corr_model.py 
