#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=calP
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/calP_%A.out
#SBATCH -t 06:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=6gb

### (1) calculate priors from ASE count tables (generated from RNA data reads)

# Load modules
module load python3
module load htslib

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
SCRIPTS=$PROJ/scripts/run_base/ase_scripts

## calculates from data
FILT=$PROJ/ase_counts_summarized

OUTPUT=$PROJ/priors_fromData
mkdir -p $OUTPUT

DESIGN=$PROJ/design_files/df_priors.csv

##### (1) calculate priors from ASE counts - data

python3 $SCRIPTS/calculate_priors_from_ase_count_tables_04amm.py \
    --output $OUTPUT \
    --design $DESIGN \
    --input $FILT \
    --debug
