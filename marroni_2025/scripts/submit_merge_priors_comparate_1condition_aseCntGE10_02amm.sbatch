#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --account=mcintyre-cid
#SBATCH --qos=mcintyre-cid
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/merge.%j.out
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=56gb



## merge priors to comparate
## merge comparates together

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
SCRIPTS=$PROJ/scripts/run_base/ase_scripts

## user must provide following design file for merging comparates:
DESIGN2=$PROJ/design_files/df_ase_compare_rils_1condition.csv

DESIGN=$PROJ/design_files/df_priors.csv
PRIORS=$PROJ/priors_fromData
FILT=$PROJ/ase_counts_summarized
#OUT=$PROJ/ase_counts_summarized_w_priors
#    mkdir -p $OUT
BAYESIN=$PROJ/bayesian_in_1condition_02amm
    mkdir -p $BAYESIN

# Load modules
module load python/3

echo "running merge priors to comparate
"

python3 ${SCRIPTS}/merge_priors_to_comparate_06amm.py \
    --output ${BAYESIN} \
    --comp ${FILT} \
    --prior ${PRIORS} \
    --design ${DESIGN}

echo "running create headers for bayesian
"
python3 ${SCRIPTS}/gen_headers_after_merge_priors_one_comparate_03amm.py \
    --output ${BAYESIN} \
    --input ${BAYESIN} \
    --design ${DESIGN2}
