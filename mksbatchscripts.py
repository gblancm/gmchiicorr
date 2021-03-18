import glob

infiles=glob.glob('./output/*corr_small.txt')

for infile in infiles:
    galaxy=infile.replace('./output/','')
    galaxy=galaxy.replace('_corr_small.txt','')
    
    with open('junk_'+galaxy+'_submit.sh', 'w') as file:
        file.write('#!/bin/bash\n')
        file.write('#\n')
        file.write('#SBATCH --job-name=gmchiicorr_'+galaxy+'\n')
        file.write('#SBATCH --nodes=1\n')
        file.write('#SBATCH --ntasks=1\n')
        file.write('#SBATCH --ntasks-per-node=24\n')
        file.write('#SBATCH --time=48:00:00\n')
        file.write('#SBATCH -p OBS,SHARED\n')
        file.write('#SBATCH --mem-per-cpu=4000\n')
        file.write('#SBATCH --output='+galaxy+'_output.out\n')
        file.write('#SBATCH --error='+galaxy+'_error.out\n')
        file.write('#SBATCH --mail-user=gblancm@carnegiescience.edu\n')
        file.write('#SBATCH --mail-type=END,FAIL\n')
        file.write('\n')
        file.write('srun python fitmodel.py '+galaxy+'\n')



