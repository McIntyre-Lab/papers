#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=id_ujc_map2one
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/idujc_map2one_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

## Description: Run ID_UJC on the self-mapped GTF

module purge

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

DESIGN_FILE=$PROJ/design_files/df_gffcomp_and_id_ujc_03ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

OUTNAME=${ARRAY[0]}
IND=${ARRAY[1]}

UJC_OUT=$PROJ/annotation_id_ujc_output
	mkdir -p $UJC_OUT

echo "Running $OUTNAME:"

# ID_UJC
python $UTIL/id_ujc_02ksb.py \
-g $IND/${OUTNAME}.gtf \
-o $UJC_OUT

echo "Num UJCs:"
NUM_UJC=$(tail -n +2 $UJC_OUT/${OUTNAME}_ujc_info.csv | wc -l | awk '{print $1}')
echo $NUM_UJC
