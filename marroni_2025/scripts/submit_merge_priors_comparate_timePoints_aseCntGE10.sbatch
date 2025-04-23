#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#BATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/merge.%j.out
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb



## merge priors to comparate
## merge comparates together

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
SCRIPTS=$PROJ/scripts/run_base/ase_scripts

## user must provide following design file for merging comparates:
DESIGN2=$PROJ/design_files/df_ase_compare_timepoints.csv

DESIGN=$PROJ/design_files/df_priors.csv
PRIORS=$PROJ/priors_fromData
FILT=$PROJ/ase_counts_summarized
BAYESIN=$PROJ/bayesian_in_01h_vs_38d
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


    python3 ${SCRIPTS}/merge_comparates_and_gen_headers_for_bayesian_w_filtering_aseGE10.py \
    --output ${BAYESIN} \
    --comp ${BAYESIN} \
    --design ${DESIGN2}
