#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=fiveSpecies_flag
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/fiveSpecies_flag_%A-%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

module purge
ml python/3.8

#export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

DESIGN_FILE=$PROJ/design_files/df_fiveSpecies_flag_02ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

GENOME=${ARRAY[0]}
SELF_MAP_LINK=(${ARRAY[1]})
GENE_KEY=${ARRAY[2]}

IND=$PROJ/mapped_ujc_id_ujc_output
OUTD=$PROJ/fiveSpecies_annotations

echo ${SELF_MAP_LINK[@]}

# Create flag file
python $UTIL/create_fiveSpecies_2_genome_flag_file_01ksb.py \
	-a ${SELF_MAP_LINK[@]} \
	-gn $GENOME \
	-i $IND \
	-gk $GENE_KEY \
	-o $OUTD/flag_fiveSpecies_2_${GENOME}_ujc.csv

