#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rlrBinXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/binXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=2-9
##SBATCH --array=2

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: Bin XLOC genes from GFFCompare output

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

if [[ $SAMPLEID == dsan* ]] ; then
        GENOME=dsan1
else
    	GENOME=dyak2
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
