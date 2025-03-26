#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pycP
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/pycP_%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=3,4,6,7,9-19,21,22,24,25,27-36,38,40,42-51,53,54,56,57,59-68

## pychopper
## rest of runs (01h_19july2022 done)

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
PYCHOP=${PROJ}/dorado_basecalling/pychopper
    mkdir -p $PYCHOP

## Get design file information (2 - 68)
DESIGN_FILE=$PROJ/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
SAMPLEID=${ARRAY[1]}
TRT=${ARRAY[4]}
TR=${ARRAY[8]}
TECH=${ARRAY[9]}
DATE=${ARRAY[10]}
BCNUM=${ARRAY[13]}

ORIG=$PROJ/dorado_basecalling/demux_1to1/${TRT}_${DATE}


    module purge
    module load pychopper/2.7.1

    echo " running pychopper $SAMPLEID
    "
    pychopper \
        -t 12 \
        -r ${PYCHOP}/report_${SAMPLEID}.pdf \
        -u ${PYCHOP}/unclass_${SAMPLEID}.fq \
        -w ${PYCHOP}/rescue_${SAMPLEID}.fq \
        ${ORIG}/SQK-PCB109_barcode${BCNUM}.fastq \
        ${PYCHOP}/fl_${SAMPLEID}.fq


    ## pychopper read counts
    echo " pychopper counts  $SAMPLEID
"

    for TYPE in fl rescue unclass; do
        INFQ=${PYCHOP}/${TYPE}_${SAMPLEID}.fq
        if [ -e "$INFQ" ]; then
            CNT=$(awk 'END {print NR/4}' ${INFQ})

            echo "sampleID,num_${TYPE}" > ${PYCHOP}/${SAMPLEID}_${TYPE}_pychop_read_cnts.csv
            echo "${SAMPLEID},${CNT}" >> ${PYCHOP}/${SAMPLEID}_${TYPE}_pychop_read_cnts.csv
        else
           echo  "File doesn't exist: $INFQ"
        fi
    done



