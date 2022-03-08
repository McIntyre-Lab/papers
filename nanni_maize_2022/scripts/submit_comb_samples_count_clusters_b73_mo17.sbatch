#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=compare_clusters_comb
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/compare_clusters_comb_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1

module purge
module load python/3.6

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/compare_b73_2_mo17/compare_clusters_in_synteny

## Combine sample clusters comparisons mapping in Mo17-B73 syntenic genes
python $PROJ/scripts/pacbio/comb_samples_count_clusters_b73_mo17.py \
    -i ${IND}
