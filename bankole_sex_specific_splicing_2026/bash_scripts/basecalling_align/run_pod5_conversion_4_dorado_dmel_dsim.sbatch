#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pod5R
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/pod5R_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=10gb
##SBATCH --partition=gpu
##SBATCH --gpus=a100:2
#SBATCH --array=9,10
date;hostname;pwd


### amm and lmm testing basecalling and demux with dorado!
## testing on patrika RILs

module load pod5/0.3.6

PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
OUT=$PROJ/dorado_basecalling
    mkdir -p $OUT

## design rows 2-10
DESIGN_FILE=$PROJ/design_files/list_dirs_w_flagDorado.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
TR=${ARRAY[1]}
TECH=${ARRAY[2]}
DATE=${ARRAY[3]}
ONT_DIR=${ARRAY[4]}

ORIG=$PROJ/original_data/${ONT_DIR}

## convert fast5 to pod5
POD5=$OUT/pod5_files_1to1
    mkdir -p $POD5

date
    cd $OUT

    pod5 convert fast5 \
        --threads 22 \
        ../original_data/${ONT_DIR}/fast5_*/*/*.fast5 \
        --output $POD5/ \
        --one-to-one ../original_data

date
