#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sqTO
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/sqTO_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=86gb
#SBATCH --array=2-9

## sqanti3QC and class count files

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
ALN=${PROJ}/aln_dorReads_2_dyak2_or_dsan1

SQANTI=${PROJ}/sqanti3QC_yak_san_dorReads_native_dyak2_dsan1_ujc
    mkdir -p $SQANTI

#CNTS=${PROJ}/gene_fsm_ism_counts_from_sqanti_classFile
#    mkdir -p $CNTS

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

ROZ=$PROJ/roz_sq_${SAMPLEID}
    mkdir -p $ROZ

#OUTD=${SQANTI}/${SAMPLEID}
#    mkdir -p ${OUTD}

if [[ $SAMPLEID == dyak* ]] ; then
    GTF=${ALN}/${SAMPLEID}_2_dyak2_noMQ0_uniq.gtf
    TAG=dyak2
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dyak21_2_dyak2_ujc.gtf
    REF_FA=/blue/mcintyre/share/references/dyak_Prin_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fa
else
    GTF=${ALN}/${SAMPLEID}_2_dsan1_noMQ0_uniq.gtf
    TAG=dsan1
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dsan11_2_dsan1_ujc.gtf
    REF_FA=/blue/mcintyre/share/references/dsan_Prin_1.1/GCF_016746245.2_Prin_Dsan_1.1_genomic.fa
fi


    module purge
    module load sqanti3/5.3.6

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID} \
        ${GTF} \
        ${REF_GTF} \
        ${REF_FA}

    ## Move all files to output directory
 #   mv ${ROZ}/*classification.txt ${OUTD}/.
 #   mv ${ROZ}/*junctions.txt ${OUTD}/.

 #   rm -r $ROZ


