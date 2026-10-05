#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=sqK
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/sqK_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=2-49

## run sqanti with self mapped serrata annotation


## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
ALN=${PROJ}/aln_ser_dorReads_2_dser1

SQANTI=${PROJ}/sqanti3QC_ser_dorReads_native_dser1_ujc
    mkdir -p $SQANTI


## Get design file information (2 - 49)
DESIGN_FILE=$PROJ/design_files/design_file_sample2BC_w_origDataPath_04amm.csv
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

ROZ=${PROJ}/rozSQ_${SAMPLEID}
    mkdir -p $ROZ

GTF=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.gtf

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dser11_2_dser1_ujc.gtf
REF_FA=/blue/mcintyre/share/references/dser1.1/GCF_002093755.2/GCF_002093755.2_Dser1.1_genomic.fa


date
    module purge
    module load sqanti3/4.2

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID} \
        ${GTF} \
        ${REF_GTF} \
        ${REF_FA}

date

## Move all files to output directory
#mv ${ROZ}/*.* ${SQANTI}/.

#rm -r $ROZ
