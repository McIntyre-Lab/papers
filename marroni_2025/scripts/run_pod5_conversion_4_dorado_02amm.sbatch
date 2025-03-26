#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pod5
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/pod5_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=40gb
##SBATCH --partition=gpu
##SBATCH --gpus=a100:2
#SBATCH --array=2
date;hostname;pwd


### amm and lmm testing basecalling and demux with dorado!
## testing on patrika RILs

module load pod5/0.3.6

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
OUT=$PROJ/dorado_basecalling
    mkdir -p $OUT

## design rows 1-6
DESIGN_FILE=$PROJ/design_files/list_ont_dirs_2_fast5.txt
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

TRT=${ARRAY[0]}
DATE=${ARRAY[1]}
ONT_DIR=${ARRAY[2]}

ORIG=$PROJ/original_data/${ONT_DIR}

## convert fast5 to pod5
POD5=$OUT/pod5_files_1to1
    mkdir -p $POD5

date
    cd $OUT

    pod5 convert fast5 \
        --threads 56 \
        --force-overwrite \
        ../original_data/${ONT_DIR}/fast5_*/*/*.fast5 \
        --output $POD5/ \
        --one-to-one ../original_data

date
