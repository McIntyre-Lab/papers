#!/bin/sh
#SBATCH --mail-user=ammmorse@ufl.edu
#SBATCH --job-name=stan2_tp
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --account=mcintyre
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/stan2_tp.%A-%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --array=2-5

# Load modules
#module purge
#module load python/3.6 htslib R/3.6

##export PATH=/blue/mcintyre/share/cegs2_MvsF_exp/conda_env/BayesASE/bin:$PATH
export PATH=/blue/mcintyre/ammorse/conda_envs/BayesAse/bin:$PATH
## 2 conditions - 01h vs 38d

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data

## Set / Create Directories and Variables
SCRIPTS=$PROJ/scripts/run_base

## directory for model testing
BAYESIN=$PROJ/bayesian_in_01h_vs_38d
BAYESOUT=$PROJ/bayesian_out_01h_vs_38d
    mkdir -p $BAYESOUT

## Set design file
DESIGN_FILE=${PROJ}/design_files/df_ase_compare_timepoints.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

COMP_1=${ARRAY[0]}
COMP_2=${ARRAY[1]}
COMPID=${ARRAY[2]}

echo "comparate 1 is $COMP_1"

## Set number of comparates to be analyzed
	## M vs F for each line = 2
COMPNUM=2

## set number of iterations = 100,000
ITER=100000

## set burn in of 10,000
WARMUP=10000

######  Run python script calling environmental bayesian model (stan2, 2 conditions)
python3 ${SCRIPTS}/NBmodel_stan2_slurm_06amm.py \
    -comparate_1 ${COMP_1} \
    -comparate_2 ${COMP_2} \
    -compID ${COMPID} \
    -datafile ${BAYESIN} \
    -datafile2 ${BAYESOUT}/ \
    -cond ${COMPNUM} \
    -workdir '/blue/mcintyre/share/BASE/scripts/stan2_scripts' \
    -routput ${BAYESOUT}/ \
    -subpath '/blue/mcintyre/share/BASE/scripts/stan2_scripts/NBModel_stan2_flex_prior_02amm.R' \
    -iterations ${ITER} \
    -warmup ${WARMUP} \
    -o ${BAYESOUT}/
