#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cmb_TRS
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/sqanti_nat_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=1-4


## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog


## Get design file information (1-4)
DESIGN_FILE=$PROJ/design_files/list_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}

OUTD=${PROJ}/sqanti_vs_native
    mkdir -p ${OUTD}


#Cat filez together keeping only header of the first file

# For classification.txt
head -n 1 $(ls $PROJ/roz_sq_${SAMPLE}*/*classification.txt | head -n 1) > ${OUTD}/${SAMPLE}_classification.txt
tail -n +2 -q $PROJ/roz_sq_${SAMPLE}*/*classification.txt >> ${OUTD}/${SAMPLE}_classification.txt

# For junctions.txt
head -n 1 $(ls $PROJ/roz_sq_${SAMPLE}*/*junctions.txt | head -n 1) > ${OUTD}/${SAMPLE}_junctions.txt
tail -n +2 -q $PROJ/roz_sq_${SAMPLE}*/*junctions.txt >> ${OUTD}/${SAMPLE}_junctions.txt



rm -r $PROJ/roz_sq_${SAMPLE}*
