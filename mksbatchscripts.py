import glob

infiles=glob.glob('./output/*corr_small.txt')

for infile in infiles:
    galaxy=infile.replace('./output/','')
    galaxy=galaxy.replace('_corr_small.txt','')
    
    with open(galaxy+'_submit.sh', 'w') as file:
        file.write('#!/bin/bash\n')
        file.write('#\n')
        file.write('#SBATCH --job-name='+galaxy+'_gmchiicorr\n')
        file.write('#SBATCH --nodes=1\n')
        file.write('#SBATCH --ntasks=1\n')
        file.write('#SBATCH --ntasks-per-node=32\n')
        file.write('#SBATCH --time=168:00:00\n')
        file.write('#SBATCH -p any\n')
        file.write('#SBATCH --mem-per-cpu=2000\n')
        file.write('#SBATCH --output='+galaxy+'_output.out\n')
        file.write('#SBATCH --error='+galaxy+'_error.out\n')
        file.write('#SBATCH --mail-user=gblancm@carnegiescience.edu\n')
        file.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
        file.write('\n')
        file.write('module load python3/3.8.5\n')
        file.write('srun python3 fitmodel.py '+galaxy+'\n')
        file.write('srun python3 plot_mcmc.py '+galaxy+'\n')

with open('submitall.sh','w') as file:
    file.write('#!/bin/bash\n')
    for infile in infiles:
        galaxy=infile.replace('./output/','')
        galaxy=galaxy.replace('_corr_small.txt','')
        file.write('sbatch '+galaxy+'_submit.sh\n')    



