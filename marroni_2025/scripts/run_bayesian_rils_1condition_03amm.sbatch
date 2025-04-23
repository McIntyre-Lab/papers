#!/bin/sh
#SBATCH --mail-user=ammmorse@ufl.edu
#SBATCH --job-name=stan2_decoupled
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --account=mcintyre-cid
#SBATCH --qos=mcintyre-cid
#SBATCH -o /blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/stan2_1cond.%A-%a.out
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=126gb
#SBATCH --array=2-9

# Load modules
#module purge
#module load python/3.6 htslib R/3.6

##export PATH=/blue/mcintyre/share/cegs2_MvsF_exp/conda_env/BayesASE/bin:$PATH
export PATH=/blue/mcintyre/ammorse/conda_envs/BayesAse/bin:$PATH

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data

## Set / Create Directories and Variables
SCRIPTS=$PROJ/scripts/run_base

BAYESIN=$PROJ/bayesian_in_1condition_02amm

BAYESOUT=$PROJ/bayesian_out_1condition_03amm
    mkdir -p $BAYESOUT

# design file with comparate and compID
DESIGN_FILE=${PROJ}/design_files/df_ase_compare_rils_1condition.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

COMP_1=${ARRAY[0]}
COMPID=${ARRAY[1]}

echo "comparate 1 is $COMP_1"
echo "compID is $COMPID"

## Set number of comparates to be analyzed
	## M vs F for each line = 2
COMPNUM=1

## set number of iterations = 100,000
ITER=100000

## set burn in of 10,000
WARMUP=10000

######  Run python script calling environmental bayesian model (stan2,  condition)
python3 ${SCRIPTS}/NBmodel_stan2_slurm_1env.py \
    -comparate_1 ${COMP_1} \
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


#/blue/mcintyre/share/BASE_mclab/galaxy/scripts/stan2_scripts
#/blue/mcintyre/share/BASE_mclab/galaxy/scripts/stan2_scripts/NBModel_stan2_flex_prior.R
