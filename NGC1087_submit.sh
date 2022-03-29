#!/bin/bash
#
#SBATCH --job-name=NGC1087_gmchiicorr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=168:00:00
#SBATCH -p any
#SBATCH --mem-per-cpu=2000
#SBATCH --output=NGC1087_output.out
#SBATCH --error=NGC1087_error.out
#SBATCH --mail-user=gblancm@carnegiescience.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load python3/3.8.5
srun python3 fitmodel.py NGC1087
srun python3 plot_mcmc.py NGC1087
