#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cat_dsim
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/cat_dsim.%j.out
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10gb

## Cat dsim ujc annotations - dsim202 + dsimmWXD
PROJ=/blue/mcintyre/share/sex_specific_splicing

IND=${PROJ}/annotation_id_ujc_output

DSIM202=${IND}/dsim202_2_dsim2_ujc.gtf

DSIMWXD=${IND}/dsimWXD_2_dsim2_ujc.gtf

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

python ${PROJ}/scripts/combine_multi_ujc_gtf_on_one_genome_01ksb.py \
    -ig ${DSIM202} ${DSIMWXD} \
    -o ${IND}/dsim202_plus_dsimWXD_2_dsim2_ujc.gtf
