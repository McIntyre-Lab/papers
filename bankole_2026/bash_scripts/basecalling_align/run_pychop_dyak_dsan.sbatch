#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pychop
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/pychop_%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
##SBATCH --array=2
#SBATCH --array=6-9

## pychopper

PROJ=/blue/mcintyre/share/transcript_ortholog
PYCHOP=${PROJ}/dorado_basecalling/pychopper
    mkdir -p $PYCHOP


## Get design file information (2 - 9, prom: 6-9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir.csv
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

ORIG=$PROJ/dorado_basecalling/demux_1to1/${TECH}_${DATE}_TR${TR}

ROZ=$PROJ/roz_dor_${SAMPLEID}
    mkdir -p $ROZ

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

    rm -r $ROZ


