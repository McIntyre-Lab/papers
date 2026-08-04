#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rmgCntXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/cntXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=2-109
##SBATCH --array=2
#SBATCH --dependency=afterok:42808173

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: count XLOC/not XLOC and output num XLOC genes and reads

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
DIR=$PROJ/gffcompare_by_TR

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_bc_design_w_origDataPath_04amm.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

if [[ $SAMPLEID == mel* ]]; then
	ALN=${PROJ}/aln_mel_dorReads_2_dmel6
	GENOME=dmel6
else
	ALN=${PROJ}/aln_sim_dorReads_2_dsim2
	GENOME=dsim2
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
