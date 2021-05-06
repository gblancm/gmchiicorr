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
        file.write('#SBATCH --ntasks-per-node=24\n')
        file.write('#SBATCH --time=12:00:00\n')
        file.write('#SBATCH -p OBS,PREEMPTION\n')
        file.write('#SBATCH --mem-per-cpu=4000\n')
        file.write('#SBATCH --output='+galaxy+'_output.out\n')
        file.write('#SBATCH --error='+galaxy+'_error.out\n')
        file.write('#SBATCH --mail-user=gblancm@carnegiescience.edu\n')
        file.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
        file.write('\n')
        file.write('srun python fitmodel.py '+galaxy+'\n')

with open('submitall.sh','w') as file:
    file.write('#!/bin/bash\n')
    for infile in infiles:
        galaxy=infile.replace('./output/','')
        galaxy=galaxy.replace('_corr_small.txt','')
        file.write('sbatch '+galaxy+'_submit.sh\n')    



