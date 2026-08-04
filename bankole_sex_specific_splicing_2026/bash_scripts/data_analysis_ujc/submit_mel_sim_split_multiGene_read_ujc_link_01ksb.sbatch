#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=ujc_lnk_R
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/ujc_lnk_R_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=126gb
#SBATCH --account=mcintyre
#SBATCH --array=1-2

module purge

## Description: split ujc link file into multigene and nonmultigene
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

OUTD=$PROJ

## Get design file information (2 rows, dmel and dsim)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=d${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

python $UTIL/output_multiGene_jxnHash_using_ujc_link_file_01ksb.py \
	-d $PROJ/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv \
	-l $PROJ/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_simpleCat.csv \
	-x ${SPECIES}_data_2_${GENOME} \
	-o $PROJ

