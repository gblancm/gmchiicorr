#!/bin/bash
#
#SBATCH --job-name=NGC1365_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=168:00:00
#SBATCH -p any
#SBATCH --mem-per-cpu=2000
#SBATCH --output=NGC1365_output.out
#SBATCH --error=NGC1365_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load python3/3.8.5
srun python3 fitmodel.py NGC1365
srun python3 plot_mcmc.py NGC1365
