#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rmgBinXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/binXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=2-109
##SBATCH --array=2

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

GTF=${DIR}/${SAMPLEID}_2_${GENOME}_updGeneID.gtf
KEY=${DIR}/${SAMPLEID}_2_${GENOME}_gffcompare_gene_key.csv

NUMREAD=$(awk '{print $12}' $GTF | sort | uniq | wc -l)
NUMGENE=$(awk '{print $10}' $GTF | sort | uniq | wc -l)

awk -F, '$3~"XLOC" {print $2}' $KEY | sort | uniq > $DIR/list_${SAMPLEID}_2_${GENOME}_XLOC_read.csv

awk '$10!~"\"XLOC_"' $GTF > ${DIR}/${SAMPLEID}_2_${GENOME}_noXLOC.gtf
awk -F, '$3!~"XLOC_"' $KEY > ${DIR}/${SAMPLEID}_2_${GENOME}_gffcompare_gene_key_noXLOC.csv

NOXLOCREADNUM=$(tail -n +2 ${DIR}/${SAMPLEID}_2_${GENOME}_gffcompare_gene_key_noXLOC.csv | awk -F, '{print $1}' | sort | uniq | wc -l)
XLOCREADNUM=$(awk '{print $1}' ${DIR}/list_${SAMPLEID}_2_${GENOME}_XLOC_read_and_gene.csv | sort | uniq | wc -l)

echo "Num reads: $NUMREAD"
echo "Num non-XLOC reads: $NOXLOCREADNUM"
echo "Num XLOC reads: $XLOCREADNUM"
SUM=$((NOXLOCREADNUM+XLOCREADNUM))
echo "Adding num non-XLOC+XLOC: $SUM"

NOXLOCGENENUM=$(tail -n +2 ${DIR}/${SAMPLEID}_2_${GENOME}_gffcompare_gene_key_noXLOC.csv | awk -F, '{print $3}' | sort | uniq | wc -l)
XLOCGENENUM=$(awk '{print $2}' ${DIR}/list_${SAMPLEID}_2_${GENOME}_XLOC_read_and_gene.csv | sort | uniq | wc -l)

echo "Num genes: $NUMGENE"
echo "Num non-XLOC genes: $NOXLOCGENENUM"
echo "Num XLOC genes: $XLOCGENENUM"
SUM=$((NOXLOCGENENUM+XLOCGENENUM))
echo "Adding num non-XLOC+XLOC: $SUM"
