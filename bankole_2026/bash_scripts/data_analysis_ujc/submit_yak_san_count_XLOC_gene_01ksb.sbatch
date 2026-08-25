#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rlrCntXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/cntXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=2-9
##SBATCH --array=2

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: count XLOC/not XLOC and output num XLOC genes and reads

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
DIR=$PROJ/gffcompare_by_TR

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
SAMPLE=${ARRAY[1]}
BC=${ARRAY[2]}
BCNUM=${ARRAY[3]}
TR=${ARRAY[4]}
TECH=${ARRAY[5]}
DATE=${ARRAY[6]}
ONT_DIR=${ARRAY[7]}

SAMPLEID=${SAMPLE}_TR${TR}

ALN=${PROJ}/aln_dorReads_2_dyak2_or_dsan1

if [[ $SAMPLEID == dsan* ]] ; then
        GENOME=dsan1
else
    	GENOME=dyak2
fi

FILE=${DIR}/${SAMPLEID}_2_${GENOME}_updGeneID.gtf
OUTFILE=$DIR/${SAMPLEID}_XLOC_cnts.csv

NUMREAD=$(awk '{print $12}' $FILE | sort | uniq | wc -l)
NUMXLOCREAD=$(awk '$10~"\"XLOC_" {print $12}' $FILE | sort | uniq | wc -l)
NUMNOXLOCREAD=$(awk '$10!~"\"XLOC_" {print $12}' $FILE | sort | uniq | wc -l)
PCTXLOCREAD=$(echo "scale=4; $NUMXLOCREAD / $NUMREAD " | bc)

NUMGENE=$(awk '{print $10}' $FILE | sort | uniq | wc -l)
NUMXLOCGENE=$(awk '$10~"\"XLOC_" {print $10}' $FILE | sort | uniq | wc -l)
NUMNOXLOCGENE=$(awk '$10!~"\"XLOC_" {print $10}' $FILE | sort | uniq | wc -l)
#PCTXLOCGENE=$(echo "scale=4; $NUMXLOCGENE / $NUMGENE " | bc)

echo "$SAMPLEID,$NUMREAD,$NUMXLOCREAD,$NUMNOXLOCREAD,$PCTXLOCREAD,$NUMGENE,$NUMXLOCGENE,$NUMNOXLOCGENE" >> $OUTFILE
