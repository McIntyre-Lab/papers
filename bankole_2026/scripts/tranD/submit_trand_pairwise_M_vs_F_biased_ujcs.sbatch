#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=trand_M_F
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/trand_M_F_%A_%a.out
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=64gb
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --array=3,6

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
IND=${PROJ}/trand_sexBias


## Get design file information
DESIGN_FILE=$PROJ/design_files/design_species_tag_directory.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
TAG=${ARRAY[1]}

OUTD=${IND}/tranD_d${SPECIES}_M_vs_F
    mkdir -p ${OUTD}

#TranD pairwise male vs female
date
trand \
${IND}/d${SPECIES}_male_ujcs_4_TranD.gtf \
${IND}/d${SPECIES}_female_ujcs_4_TranD.gtf \
-o $OUTD \
-1 d${SPECIES}_male \
-2 d${SPECIES}_female \
-f \
-n $SLURM_CPUS_ON_NODE
date
