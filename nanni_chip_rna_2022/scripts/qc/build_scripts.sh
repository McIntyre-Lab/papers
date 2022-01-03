#!/bin/sh

### Variables from master submission script:
# ${EMAIL}        your email address
# ${PROJ}         project folder
# ${ORIG}         original data folder
# ${FASTQ}        Un-gzipped FASTQ folder
# ${QC}           QC folder
# ${LOGDIR}       PBS logs folder
# ${JOBARRAYIDS}  job array (e.g. 2-97)
# ${QCPIPE}       Folder of the QC pipeline templates and scripts
# ${PE}           Flag if including PE alignments in pipeline


######## GENERATE QSUBS
SCRIPTS=$PROJ/scripts/qc
   if [ ! -e $SCRIPTS ]; then mkdir -p $SCRIPTS; fi
TEMPLATES=$QCPIPE/templates

### Un-Gzip FASTQ files qsub

JOB0_NAME=gunzip_fastq
JOB0_TIME=1:00:00
JOB0_NODES=1
JOB0_MEM=1gb
JOB0_SCRIPT=$SCRIPTS/gunzip_fastq.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB0_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB0_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB0_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB0_MEM}/g" > ${JOB0_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB0_SCRIPT}
echo "#" >> ${JOB0_SCRIPT}
echo >> ${JOB0_SCRIPT}
echo "## Set directories" >> ${JOB0_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB0_SCRIPT}
echo "FASTQ=${FASTQ}" >> ${JOB0_SCRIPT}
echo "if [ ! -e \$FASTQ ]; then mkdir -p \$FASTQ; fi" >> ${JOB0_SCRIPT}
echo "ORIG=${ORIG}" >> ${JOB0_SCRIPT}
echo >> ${JOB0_SCRIPT}
echo >> ${JOB0_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB0_SCRIPT}
echo >> ${JOB0_SCRIPT}
echo >> ${JOB0_SCRIPT}
echo "## Gunzip files" >> ${JOB0_SCRIPT}
echo "gunzip -c \$ORIG/\${LOCATION}/\${READ}.gz > \$FASTQ/\${READ}" >> ${JOB0_SCRIPT}


source $QCPIPE/scripts/build_scripts/build_fastqc_scripts.sh
source $QCPIPE/scripts/build_scripts/build_ercc_scripts.sh
#source $QCPIPE/scripts/build_scripts/build_cleanup_script.sh

