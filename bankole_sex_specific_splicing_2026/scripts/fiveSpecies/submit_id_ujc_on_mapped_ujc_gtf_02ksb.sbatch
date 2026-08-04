#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=id_ujc_map2each
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/idujc_map2each_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=1-24
#SBATCH --qos=mcintyre-b

## Description: ID_UJC for the mapped UJC GTF (ujc gtfs mapped to other genomes)

module purge

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

# list of mapped UJC GTFs
DESIGN_FILE=$PROJ/design_files/df_noHeader_id_ujc_mapped_ujc_gtf_02ksb.csv

OUTD=$PROJ/mapped_ujc_id_ujc_output
	mkdir -p $OUTD

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

INFILE=${ARRAY[0]}

echo "File: $INFILE"

# ID_UJC
python $UTIL/id_ujc_02ksb.py \
-g $INFILE \
-o $OUTD
