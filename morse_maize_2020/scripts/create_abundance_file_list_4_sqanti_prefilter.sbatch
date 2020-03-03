#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=abun
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/abun_%A.out
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=1


### create a file containing the abundance.txt files for each sample 
### input into sqanti prefiltering 



### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

cd $PROJ/tofu2_b73

find -maxdepth 2 -type f -name "*abundance*" | sed 's@./@${PROJ}/tofu2_b73/@' | paste -sd, >$PROJ/tofu2_b73/abundance_file_list.csv
