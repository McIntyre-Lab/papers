#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=compare_clusters
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/compare_clusters_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1,5,11

## Run for just mo17 samples (array 1,5,11)

module purge
module load python/3.6

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
SQANTIb=$PROJ/compare_b73_2_mo17/sqanti_b73
TOFUb=$PROJ/compare_b73_2_mo17/tofu2_b73
SQANTIm=$PROJ/compare_b73_2_mo17/sqanti_mo17_cau
TOFUm=$PROJ/compare_b73_2_mo17/tofu2_mo17_cau
OUTD=$PROJ/compare_b73_2_mo17/compare_clusters_in_synteny
    mkdir -p ${OUTD}


## Get B73-Mo17 CAU COGE synteny file
SYN=/blue/mcintyre/share/references/maize_mo17/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv

## Get supp. file of 12604 genes of interest
SUPP=$PROJ/supp_from_resubmission_2021/Supplementary_File_5.csv

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

## Set input directories to the correct species
SbDIR=${SQANTIb}/${GENO}/${TRT}
TbDIR=${TOFUb}/${SAMPLE_NAME}
SmDIR=${SQANTIm}/${GENO}/${TRT}
TmDIR=${TOFUm}/${SAMPLE_NAME}

## Compare clusters mapping in Mo17-B73 syntenic genes
python $PROJ/scripts/pacbio/compare_clusters_b73_mo17_cau_syntenic.py \
    -i ${SUPP} \
    --synteny ${SYN} \
    -tB73 ${TbDIR}/${SAMPLE_NAME}.collapsed.read_stat.txt \
    -tMo17 ${TmDIR}/${SAMPLE_NAME}.collapsed.read_stat.txt \
    -sB73 ${SbDIR}/${SAMPLE_NAME}_classification.txt \
    -sMo17 ${SmDIR}/${SAMPLE_NAME}_classification.txt \
    -p ${OUTD}/${SAMPLE_NAME}_b73_mo17CAU_cluster_compare
