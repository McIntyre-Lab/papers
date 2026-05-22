#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=flgFileIR
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/flgFileIR_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=0-4

module purge

export PATH=/ufgi-grimshawi/TB14/k.bankole/conda_env/TranD/bin:$PATH
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/

GENOME_ARR=(dmel6 dsim2 dsan1 dyak2 dser1)

GENOME=${GENOME_ARR[$SLURM_ARRAY_TASK_ID]}

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

IND=$PROJ/fiveSpecies_annotations

python $UTIL/add_wIR_flag_to_fiveSpecies_flag_file_01ksb.py \
	-l $IND/list_noHeader_fiveSpecies_2_${GENOME}_genes_missing_from_event_analysis.csv \
	-f $IND/flag_fiveSpecies_2_${GENOME}_ujc.csv \
	-o $IND/flag_fiveSpecies_2_${GENOME}_ujc_w_IR_flag.csv

