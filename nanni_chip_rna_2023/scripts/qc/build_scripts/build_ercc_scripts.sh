#!/bin/bash


###### ERCC alignments (single ended)

JOB9_NAME=ercc_aln
JOB9_TIME=2:00:00
JOB9_NODES=4
JOB9_MEM=6gb
JOB9_SCRIPT=$SCRIPTS/ercc_aln_se.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB9_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB9_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB9_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB9_MEM}/g" > ${JOB9_SCRIPT}
echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB9_SCRIPT}
echo "#" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo 'mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' >> ${JOB9_SCRIPT}
echo 'export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "# Load modules" >> ${JOB9_SCRIPT}
echo "module load bowtie/0.12.9" >> ${JOB9_SCRIPT}
echo "module load python/2.7.6" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Store number of processors that will be used, i.e. the number of files we will split into to run LAST on" >> ${JOB9_SCRIPT}
echo "    NUMPROCS=4" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Set Directories" >> ${JOB9_SCRIPT}
echo "      PROJ=${PROJ}" >> ${JOB9_SCRIPT}
echo "      ORIG=${FASTQ}" >> ${JOB9_SCRIPT}
echo "      QC=${QC}" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Set Different References" >> ${JOB9_SCRIPT}
echo "    REF=/ufrc/mcintyre/share/references/ERCC_Ambion/ERCC92" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Create Necessary Folders" >> ${JOB9_SCRIPT}
echo "    OUTPUT=\$QC/ercc/aln_ERCC_all_se" >> ${JOB9_SCRIPT}
echo "    if [ ! -e \$OUTPUT ] ; then mkdir -p \$OUTPUT; fi" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "    # Create JOB LOG directory and start log" >> ${JOB9_SCRIPT}
echo "        LOGS=\$OUTPUT/job_logs" >> ${JOB9_SCRIPT}
echo "        ALN_ERROR_LOG=\$LOGS/size_errors.txt" >> ${JOB9_SCRIPT}
echo "        if [ ! -d \$LOGS ]; then mkdir -p \$LOGS; fi" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "        MYLOG=\$LOGS/\${NAME}.log" >> ${JOB9_SCRIPT}
echo '        printf "`date` $NAME PBS_ARRAYID:$PBS_ARRAYID HOSTNAME:$HOSTNAME\n" > $MYLOG' >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "    # Create ALN LOG directory" >> ${JOB9_SCRIPT}
echo "        ALNLOGS=\$OUTPUT/aln_logs" >> ${JOB9_SCRIPT}
echo "        if [ ! -d \$ALNLOGS ]; then mkdir -p \$ALNLOGS; fi" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "    # Create UNALN READ directory" >> ${JOB9_SCRIPT}
echo "        UNALNDIR=\$OUTPUT/unaln_reads" >> ${JOB9_SCRIPT}
echo "        if [ ! -d \$UNALNDIR ]; then mkdir -p \$UNALNDIR; fi" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "### FUNCTIONS FOR ALIGNMENT PIPELINE ###" >> ${JOB9_SCRIPT}
echo "# I have created a separate file that hold the alignment functions" >> ${JOB9_SCRIPT}
echo "source $QCPIPE/scripts/alignment_functions.sh" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Start Alignment Pipeline" >> ${JOB9_SCRIPT}
echo '    printf "<------------------- STARTING SE alignment process for $READ1 [`date`] ------------------->\n" >> $MYLOG' >> ${JOB9_SCRIPT}
echo "    READS=\$ORIG/\${READ}" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "    qual=\`python /ufrc/mcintyre/share/python.git/identify_quality.py -i \$READS\`" >> ${JOB9_SCRIPT}
echo '    if [ $qual == "phred64" ]; ' >> ${JOB9_SCRIPT}
echo "    then " >> ${JOB9_SCRIPT}
echo "        # set to old illumina quality scores phred64/solexa 1.3" >> ${JOB9_SCRIPT}
echo '        btqual="--phred64-quals"' >> ${JOB9_SCRIPT}
echo '        lastqual="3"' >> ${JOB9_SCRIPT}
echo "    else" >> ${JOB9_SCRIPT}
echo "        # change to sanger format which is what all new illumina data is" >> ${JOB9_SCRIPT}
echo '        btqual="--phred33-quals"' >> ${JOB9_SCRIPT}
echo '        lastqual="1"' >> ${JOB9_SCRIPT}
echo "    fi" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "    bowtie_se_uniq" >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Combine all Sam files" >> ${JOB9_SCRIPT} 
echo '    echo "START Combine SAM files">>$MYLOG' >> ${JOB9_SCRIPT}
echo "    cat *.sam >\$OUTPUT/\${NAME}.sam 2>>\$MYLOG" >> ${JOB9_SCRIPT}
echo '    echo "FINISH Combining SAM files [`date`]" >>$MYLOG' >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo "## Combine all Unaln FASTQ " >> ${JOB9_SCRIPT}
echo '    echo "START Combine Unaln FQ">>$MYLOG' >> ${JOB9_SCRIPT}
echo "    cat *_unaln_bt.fq >\$UNALNDIR/\${NAME}_unaln.fq 2>>\$MYLOG" >> ${JOB9_SCRIPT}
echo '    echo "FINISH Combining SAM files [`date`]" >>$MYLOG' >> ${JOB9_SCRIPT}
echo >> ${JOB9_SCRIPT}
echo 'echo "Script complete, [`date`]" >>$MYLOG' >> ${JOB9_SCRIPT}
echo "### Clean up temp dir" >> ${JOB9_SCRIPT}
echo "cd \$PROJ/scripts/SLURM_LOGS" >> ${JOB9_SCRIPT}
echo 'rm -f $TMPDIR\*' >> ${JOB9_SCRIPT}

###### ERCC mpileups (single ended)

JOB10_NAME=ercc_mpileups
JOB10_TIME=2:00:00
JOB10_NODES=1
JOB10_MEM=6gb
JOB10_SCRIPT=$SCRIPTS/ercc_mpileups_se.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB10_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB10_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB10_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB10_MEM}/g" > ${JOB10_SCRIPT}
echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB10_SCRIPT}
echo "#" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "# Load modules" >> ${JOB10_SCRIPT}
echo "module load samtools/1.3.1 " >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "## Set Directories" >> ${JOB10_SCRIPT}
echo "      PROJ=${PROJ}" >> ${JOB10_SCRIPT}
echo "      QC=${QC}" >> ${JOB10_SCRIPT}
echo "      INPUT=\$QC/ercc/aln_ERCC_all_se" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "# Reference" >> ${JOB10_SCRIPT}
echo "REF=/ufrc/mcintyre/share/references/ERCC_Ambion/ERCC92.fa" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "SAM=\$INPUT/\${NAME}.sam" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "#### create output directory" >> ${JOB10_SCRIPT}
echo "OUTPUT=\$QC/ercc/ercc_se_mpileups" >> ${JOB10_SCRIPT}
echo "    if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB10_SCRIPT}
echo "LOGS=\$OUTPUT/logs" >> ${JOB10_SCRIPT}
echo "    if [ ! -e \$LOGS ]; then mkdir \$LOGS; fi" >> ${JOB10_SCRIPT}
echo "MYLOG=\$LOGS/\${NAME}_mpileup.log" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "BAMOUT=\$QC/ercc/bam_files_ercc" >> ${JOB10_SCRIPT}
echo "    if [ ! -e \$BAMOUT ]; then mkdir -p \$BAMOUT; fi" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo "#### Convert SAM to BAM and make mpileups" >> ${JOB10_SCRIPT}
echo "    BAM=\$BAMOUT/\${NAME}" >> ${JOB10_SCRIPT}
echo >> ${JOB10_SCRIPT}
echo '    printf "<-------------------- Convert SAM to BAM -------------------->" >> "${MYLOG}"' >> ${JOB10_SCRIPT}
echo '    echo `date`": Starting SAM to BAM conversion" >> "${MYLOG}"' >> ${JOB10_SCRIPT}
echo '    samtools view -ut $REF.fai -o $BAM.bam $SAM 2>> "${MYLOG}" ' >> ${JOB10_SCRIPT} 
echo '    samtools sort -T ${BAM}.tmp.sorted -o $BAM.sorted.bam $BAM.bam 2>> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo '    samtools index $BAM.sorted.bam >> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo >> ${JOB10_SCRIPT}
echo '    echo `date`": Finished SAM to BAM conversion" >> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo >> ${JOB10_SCRIPT}
echo '#### Make mpielup' >> ${JOB10_SCRIPT} 
echo '    PILEUP=$OUTPUT/${NAME}.mpileup' >> ${JOB10_SCRIPT} 
echo >> ${JOB10_SCRIPT}
echo '    printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo '    echo `date`": Generating pileup" >> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo '    samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> "${MYLOG}"' >> ${JOB10_SCRIPT} 
echo >> ${JOB10_SCRIPT}
echo 'echo `date`": Finished Script complete" >> "${MYLOG}"' >> ${JOB10_SCRIPT} 

###### ERCC coverage counts (single ended)

JOB11_NAME=ercc_coverage
JOB11_TIME=0:30:00
JOB11_NODES=1
JOB11_MEM=20gb
JOB11_SCRIPT=$SCRIPTS/ercc_coverage_counts.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB11_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB11_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB11_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB11_MEM}/g" > ${JOB11_SCRIPT}
echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB11_SCRIPT}
echo "#" >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo 'mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' >> ${JOB11_SCRIPT}
echo 'export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "# Load modules" >> ${JOB11_SCRIPT}
echo "module load python " >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "## Set Directories" >> ${JOB11_SCRIPT}
echo "      PROJ=${PROJ}" >> ${JOB11_SCRIPT}
echo "      QC=${QC}" >> ${JOB11_SCRIPT}
echo "      RPKM=$QCPIPE/scripts/rpkm_calculate.py" >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "#### Make Output Directory" >> ${JOB11_SCRIPT}
echo "        OUTPUT=\$QC/ercc/coverage_counts_ercc" >> ${JOB11_SCRIPT}
echo "        if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "    # Create LOG directory and start log" >> ${JOB11_SCRIPT}
echo "        LOGS=\$OUTPUT/logs" >> ${JOB11_SCRIPT}
echo "        if [ ! -e \$LOGS ]; then mkdir -p \$LOGS; fi" >> ${JOB11_SCRIPT}
echo "        MYLOG=\$LOGS/\${FILE}.log" >> ${JOB11_SCRIPT}
echo '        printf "`date` $FILE PBS_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}" ' >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "#### COVERAGE COUNTS" >> ${JOB11_SCRIPT}
echo "    BED=/ufrc/mcintyre/share/references/ERCC_Ambion/ERCC92_v2.bed" >> ${JOB11_SCRIPT}
echo "    SAM=\$QC/ercc/aln_ERCC_all_se/\${NAME}.sam" >> ${JOB11_SCRIPT}
echo "    MPILEUP=\$QC/ercc/ercc_se_mpileups/\${NAME}.mpileup" >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo '    echo "Starting Coverage Counts for $NAME `date`" > "${MYLOG}"' >> ${JOB11_SCRIPT}
echo '    python $RPKM -b $BED -m ${MPILEUP} -s ${SAM} -n ${NAME} --cv -g "${MYLOG}" -o $OUTPUT/cvrg_cnts_${NAME}.csv' >> ${JOB11_SCRIPT}
echo '    echo "Finished Coverage Counts for $NAME `date`" >> "${MYLOG}"' >> ${JOB11_SCRIPT}
echo >> ${JOB11_SCRIPT}
echo "### Clean up temp dir" >> ${JOB11_SCRIPT}
echo "cd \$PROJ/scripts/SLURM_LOGS" >> ${JOB11_SCRIPT}
echo 'rm -f $TMPDIR\*' >> ${JOB11_SCRIPT}

###### ERCC concentration plots

JOB12_NAME=ercc_conc_plots
JOB12_TIME=2:00:00
JOB12_NODES=1
JOB12_MEM=20gb
JOB12_SCRIPT=$SCRIPTS/ercc_conc_plots.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB12_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB12_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB12_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB12_MEM}/g" > ${JOB12_SCRIPT}
echo "#SBATCH --array=1" >> ${JOB12_SCRIPT}
echo "#" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "# Load modules" >> ${JOB12_SCRIPT}
echo "module load sas/9.4 " >> ${JOB12_SCRIPT}
echo "module load R " >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "## Set Directories" >> ${JOB12_SCRIPT}
echo "      PROJ=${PROJ}" >> ${JOB12_SCRIPT}
echo "      QC=${QC}" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "      SASPROG=\$PROJ/scripts/qc/ercc_adj_conc.sas" >> ${JOB12_SCRIPT}
echo "      SASLIB=qc" >> ${JOB12_SCRIPT}
echo '      QCSASPATH=$( echo ${QC} | sed -e "s/\//\\\\\//g")' >> ${JOB12_SCRIPT}
echo '      ERCCCVGPATH=$( echo ${QC}/ercc/coverage_counts_ercc | sed -e "s/\//\\\\\//g")' >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "      CONCPLOT=$QCPIPE/scripts/concentration_plots_ERCC.R" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "      ERCC_CVG=\$QC/ercc/ercc_coverage_counts.csv" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "## Build SAS program" >> ${JOB12_SCRIPT}
echo "less $TEMPLATES/ercc_adj_conc.sas | " >> ${JOB12_SCRIPT}
echo '    sed -e "s/_SASLIB_/${SASLIB}/g" | ' >> ${JOB12_SCRIPT}
echo '    sed -e "s/_PROJPATH_/${QCSASPATH}/g" | ' >> ${JOB12_SCRIPT}
echo '    sed -e "s/_CVGPATH_/${ERCCCVGPATH}/g" > ${SASPROG} ' >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "## Combine ERCC coverage counts" >> ${JOB12_SCRIPT}
echo 'cd $QC/ercc/coverage_counts_ercc' >> ${JOB12_SCRIPT}
echo 'FLAG=0' >> ${JOB12_SCRIPT}
echo 'for FILE in *.csv' >> ${JOB12_SCRIPT}
echo 'do' >> ${JOB12_SCRIPT}
echo '    if [ $FLAG == 0 ]' >> ${JOB12_SCRIPT}
echo '    then' >> ${JOB12_SCRIPT}
echo '    cat $FILE > ${ERCC_CVG}' >> ${JOB12_SCRIPT}
echo '    FLAG=1' >> ${JOB12_SCRIPT}
echo '    else' >> ${JOB12_SCRIPT}
echo '    tail -n +2 $FILE >> ${ERCC_CVG}' >> ${JOB12_SCRIPT}
echo '    fi' >> ${JOB12_SCRIPT}
echo 'done' >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "## Run SAS" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "LOGS=\$QC/ercc/sas_logs" >> ${JOB12_SCRIPT}
echo "if [ ! -e \$LOGS ]; then mkdir -p \$LOGS; fi" >> ${JOB12_SCRIPT}
echo "SASDATA=\$QC/ercc/sas_data" >> ${JOB12_SCRIPT}
echo "if [ ! -e \$SASDATA ]; then mkdir -p \$SASDATA; fi" >> ${JOB12_SCRIPT}
echo "TMPDIR=\$QC/sas_temp" >> ${JOB12_SCRIPT}
echo "if [ ! -e \$TMPDIR ]; then mkdir -p \$TMPDIR; fi" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "sas -log \${SASPROG}.log -work \$TMPDIR -sysin \$SASPROG -sysparm \${ERCC_CVG}" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "rm \$TMPDIR/*" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "## Make concentration plots" >> ${JOB12_SCRIPT}
echo >> ${JOB12_SCRIPT}
echo "Rscript \$CONCPLOT \$QC/ercc/ercc_cvg_cnts_adj_conc.csv \$QC/ercc/ERCC_conc_plots.pdf" >> ${JOB12_SCRIPT}

