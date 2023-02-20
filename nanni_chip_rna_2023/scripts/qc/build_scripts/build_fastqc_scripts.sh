#!/bin/sh

### Build FASTQC qsub
JOB1_NAME=fastqc
JOB1_TIME=2:00:00
JOB1_NODES=1
JOB1_MEM=900mb
JOB1_SCRIPT=$SCRIPTS/run_fastqc.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB1_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB1_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB1_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB1_MEM}/g" > ${JOB1_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB1_SCRIPT}
echo "#" >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo "# Load modules" >> ${JOB1_SCRIPT}
echo "module load fastqc/0.10.1" >> ${JOB1_SCRIPT}
echo "module load htmldoc" >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo "#Set directories" >> ${JOB1_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB1_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB1_SCRIPT}
echo "QC=${QC}" >> ${JOB1_SCRIPT}
echo "OUTPUT=\$QC/fastqc/fastqc_files" >> ${JOB1_SCRIPT}
echo "if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB1_SCRIPT}
echo "PDF=\$QC/fastqc/fastqc_pdfs" >> ${JOB1_SCRIPT}
echo "if [ ! -e \$PDF ]; then mkdir -p \$PDF; fi" >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo "#Run FASTQC" >> ${JOB1_SCRIPT}
echo "fastqc \${ORIG}/\${READ} --outdir=\$OUTPUT" >> ${JOB1_SCRIPT}
echo >> ${JOB1_SCRIPT}
echo "#Convert to PDFs" >> ${JOB1_SCRIPT}
echo "READFILE=\$(basename \${READ} .fastq)" >>  ${JOB1_SCRIPT}
echo "htmldoc --webpage --fontsize 7 --browserwidth 800 --header . \$OUTPUT/\${READFILE}_fastqc/fastqc_report.html -f \$PDF/\${NAME}.pdf" >> ${JOB1_SCRIPT}


### Build FASTQC summary qsub

JOB2_NAME=fastqc_summary
JOB2_TIME=2:00:00
JOB2_NODES=1
JOB2_MEM=1gb
JOB2_SCRIPT=$SCRIPTS/run_fastqc_summary.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB2_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB2_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB2_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB2_MEM}/g" > ${JOB2_SCRIPT}
echo "#SBATCH --array=1" >> ${JOB2_SCRIPT}
echo "#" >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo "#Set directories" >> ${JOB2_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB2_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB2_SCRIPT}
echo "QC=${QC}" >> ${JOB2_SCRIPT}
echo "DATAIN=\$QC/fastqc/fastqc_files" >> ${JOB2_SCRIPT}
echo "DATAOUT=\$QC/fastqc/fastqc_summary" >> ${JOB2_SCRIPT}
echo "if [ ! -e \$DATAOUT ]; then mkdir -p \$DATAOUT; fi" >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo "STATSOUT=\$DATAOUT/files" >> ${JOB2_SCRIPT}
echo "if [ ! -e \$STATSOUT ]; then mkdir -p \$STATSOUT; fi" >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo "PASSFAIL=$QCPIPE/scripts/fastqc_pass-fail_list_rlb.pl" >> ${JOB2_SCRIPT}
echo "STATS=$QCPIPE/scripts/fastqc_stats_v2_jmf.pl" >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo >> ${JOB2_SCRIPT}
echo 'echo "Sample_id,basic_stats,base_qual,seq_qual,base_content,base_GC,seq_GC,base_N,seq_length,seq_dup,seq_overrep,kmer,percent_pass" > $DATAOUT/fastqc_pass-fail_summary.csv' >> ${JOB2_SCRIPT}
echo "cd \$DATAIN" >> ${JOB2_SCRIPT}
echo  >> ${JOB2_SCRIPT}
echo "for SAMPLE in \$( find ./ -maxdepth 1 -type d -print | cut -d'/' -f2)" >> ${JOB2_SCRIPT}
echo "    do" >> ${JOB2_SCRIPT}
echo "    FILE=\$DATAIN/\$SAMPLE/summary.txt" >> ${JOB2_SCRIPT}
echo "    perl \$FASTQC \$SAMPLE \$FILE  >> \$DATAOUT/fastqc_pass-fail_summary.csv" >> ${JOB2_SCRIPT}
echo "done" >> ${JOB2_SCRIPT}
echo  >> ${JOB2_SCRIPT}
echo "for DIR in \$( find ./ -maxdepth 1 -type d | cut -f2 -d'/')" >> ${JOB2_SCRIPT}
echo "    do" >> ${JOB2_SCRIPT}
echo "    cd \$DATAIN/\$DIR" >> ${JOB2_SCRIPT}
echo "    NAME=\`basename \$DIR _fastqc\`" >> ${JOB2_SCRIPT}
echo "    perl \$STATS \$NAME \$STATSOUT" >> ${JOB2_SCRIPT}
echo "done" >> ${JOB2_SCRIPT}


### Build identify duplicates qsub

JOB3_NAME=dupcounts
JOB3_TIME=2:00:00
JOB3_NODES=1
JOB3_MEM=20gb
JOB3_SCRIPT=$SCRIPTS/run_dup_counts.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB3_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB3_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB3_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB3_MEM}/g" > ${JOB3_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB3_SCRIPT}
echo "#" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "# Load modules" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "module load python/2.7.6" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "#Set directories" >> ${JOB3_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB3_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB3_SCRIPT}
echo "QC=${QC}" >> ${JOB3_SCRIPT}
echo "MCPYTHON=/ufrc/mcintyre/share/python.git" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "OUTPUT=\$QC/fastq_split_dups" >> ${JOB3_SCRIPT}
echo "if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "OUTCOUNTS=\$QC/duplicate_counts" >> ${JOB3_SCRIPT}
echo "if [ ! -e \$OUTCOUNTS ]; then mkdir -p \$OUTCOUNTS; fi" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "OUTLOGS=\$QC/duplicate_counts/logs" >> ${JOB3_SCRIPT}
echo "if [ ! -e \$OUTLOGS ]; then mkdir -p \$OUTLOGS; fi" >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo >> ${JOB3_SCRIPT}
echo "\$MCPYTHON/fastqSplitDups.py -r1 \$ORIG/\${READ} --outdir \$OUTPUT -o \$OUTCOUNTS/\${NAME}_duplicate_counts_summary.csv -t \$OUTCOUNTS/\${NAME}_001_duplicate_sequences.tsv -g \$OUTLOGS/\${NAME}_001_fastqSplitDups.log" >> ${JOB3_SCRIPT}

### Build identify homopolymers qsub

JOB4_NAME=homopolymers
JOB4_TIME=1:00:00
JOB4_NODES=1
JOB4_MEM=6gb
JOB4_SCRIPT=$SCRIPTS/run_identify_homopolymers.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB4_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB4_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB4_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB4_MEM}/g" > ${JOB4_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB4_SCRIPT}
echo "#" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "# Load modules" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "module load python/2.7.6" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "#Set directories" >> ${JOB4_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB4_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB4_SCRIPT}
echo "QC=${QC}" >> ${JOB4_SCRIPT}
echo "MCPYTHON=/ufrc/mcintyre/share/python.git" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "OUTPUT=\$QC/homopolymer_files" >> ${JOB4_SCRIPT}
echo "if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "OUTPNG=\$QC/homopolymer_pngs" >> ${JOB4_SCRIPT}
echo "if [ ! -e \$OUTPNG ]; then mkdir -p \$OUTPNG; fi" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "OUTLOG=\$QC/homopolymer_logs" >> ${JOB4_SCRIPT}
echo "if [ ! -e \$OUTLOG ]; then mkdir -p \$OUTLOG; fi" >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo >> ${JOB4_SCRIPT}
echo "cd \${ORIG}" >> ${JOB4_SCRIPT}
echo "# run script" >> ${JOB4_SCRIPT}
echo "        LOG=\$OUTLOG/\${NAME}.log" >> ${JOB4_SCRIPT}
echo "        PNG=\$OUTPNG/\${NAME}.png" >> ${JOB4_SCRIPT}
echo "        python \$MCPYTHON/identify_homopolymers.py --input \${ORIG}/\${READ} --fastq -a --out \$OUTPUT/\${NAME}.csv --log \$LOG --plot \$PNG" >> ${JOB4_SCRIPT}

### Build fastq2fasta qsub

JOB5_NAME=fastq2fasta
JOB5_TIME=2:00:00
JOB5_NODES=1
JOB5_MEM=12gb
JOB5_SCRIPT=$SCRIPTS/run_fastq2fasta.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB5_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB5_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB5_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB5_MEM}/g" > ${JOB5_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB5_SCRIPT}
echo "#" >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo "# Load modules" >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo "module load python/2.7.6" >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo "#Set directories" >> ${JOB5_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB5_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB5_SCRIPT}
echo "QC=${QC}" >> ${JOB5_SCRIPT}
echo "MCPYTHON=/ufrc/mcintyre/share/python.git" >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo "OUTPUT=\$QC/fasta_reads" >> ${JOB5_SCRIPT}
echo "if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo >> ${JOB5_SCRIPT}
echo "# run script " >> ${JOB5_SCRIPT}
echo "        python \$MCPYTHON/fastq2fasta.py -i \${ORIG}/\${READ} -o \$OUTPUT/\${READ}.fa" >> ${JOB5_SCRIPT}


### Build blat adapters qsub

JOB6_NAME=illumina_blat
JOB6_TIME=6:00:00
JOB6_NODES=1
JOB6_MEM=8gb
JOB6_SCRIPT=$SCRIPTS/run_blat_illumina_qsub.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB6_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB6_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB6_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB6_MEM}/g" > ${JOB6_SCRIPT}

echo "#SBATCH --array=${JOBARRAYIDS}" >> ${JOB6_SCRIPT}
echo "#" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "mkdir -p tmp/\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}" >> ${JOB6_SCRIPT}
echo 'export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}' >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "# Load modules" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "module load blat" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}

echo "#Set directories" >> ${JOB6_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB6_SCRIPT}
echo "QC=${QC}" >> ${JOB6_SCRIPT}
echo "ORIG=\${QC}/fasta_reads" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "OUTDIR=\$QC/qc_adapters" >> ${JOB6_SCRIPT}
echo "if [ ! -e \$OUTDIR ]; then mkdir -p \$OUTDIR; fi" >> ${JOB6_SCRIPT}
echo "OUTFILES=\$OUTDIR/files" >> ${JOB6_SCRIPT}
echo "if [ ! -e \$OUTFILES ]; then mkdir -p \$OUTFILES; fi" >> ${JOB6_SCRIPT}
echo "OUTLOGS=\$OUTDIR/logs" >> ${JOB6_SCRIPT}
echo "if [ ! -e \$OUTLOGS ]; then mkdir -p \$OUTLOGS; fi" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
less $TEMPLATES/design_array_slurm.txt >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "## Initialize Log" >> ${JOB6_SCRIPT}
echo "    MYLOG=\$OUTLOGS/\${NAME}_adapter.log" >> ${JOB6_SCRIPT}
echo '    printf "`date` $NAME ARRAY_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG' >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo "## Set Reference" >> ${JOB6_SCRIPT}
echo "    REF=/ufrc/mcintyre/share/references/illumina_adapters/illumina_adapter_list.fasta" >> ${JOB6_SCRIPT}
echo '    echo "file,total_num_reads,num_reads_w_adapter,per_adapter" > $OUTFILES/${READ}.csv' >> ${JOB6_SCRIPT}
echo "## Start BLAT" >> ${JOB6_SCRIPT}
echo "   cd \$ORIG" >> ${JOB6_SCRIPT}
echo "   FILE=\${READ}.fa" >> ${JOB6_SCRIPT}
echo "   if [ -s \$FILE ]" >> ${JOB6_SCRIPT}
echo "  then" >> ${JOB6_SCRIPT}
echo "      SAMPLE=\`basename \$FILE .fa\`" >> ${JOB6_SCRIPT}
echo '      printf "`date` Starting BLAT $FILE\n" >> $MYLOG' >> ${JOB6_SCRIPT}
echo "      blat \$REF \$FILE  \$TMPDIR/\${READ}.psl 2>> \$MYLOG" >> ${JOB6_SCRIPT}
echo '      printf "`date` Finished BLAT $FILE\n" >> $MYLOG' >> ${JOB6_SCRIPT}
echo '      TOTAL=`grep ">" $FILE | wc -l`' >> ${JOB6_SCRIPT}
echo "      COUNT=\`cut -f 10 \$TMPDIR/\${READ}.psl | sort | uniq | wc -l\`" >> ${JOB6_SCRIPT}
echo '      PER=`echo "scale=4; $COUNT / $TOTAL * 100" | bc`' >> ${JOB6_SCRIPT}
echo '      echo "$SAMPLE,$TOTAL,$COUNT,$PER" >> $OUTFILES/${READ}.csv' >> ${JOB6_SCRIPT}
echo "  else" >> ${JOB6_SCRIPT}
echo '      echo "`date` File was empty." >> $MYLOG' >> ${JOB6_SCRIPT}
echo "  fi" >> ${JOB6_SCRIPT}
echo >> ${JOB6_SCRIPT}
echo 'printf "`date` Script complete" >> $MYLOG' >> ${JOB6_SCRIPT}
echo "### Clean up temp dir" >> ${JOB6_SCRIPT}
echo "cd \$PROJ/scripts/SLURM_LOGS" >> ${JOB6_SCRIPT}
echo 'rm -f $TMPDIR\*' >> ${JOB6_SCRIPT}

### Summarize QC files

JOB7_NAME=summarize_qc
JOB7_TIME=1:00:00
JOB7_NODES=1
JOB7_MEM=4gb
JOB7_SCRIPT=$SCRIPTS/summarize_qc_outputs.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB7_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB7_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB7_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB7_MEM}/g" > ${JOB7_SCRIPT}
echo "#SBATCH --array=1" >> ${JOB7_SCRIPT}
echo "#" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "#Set directories" >> ${JOB7_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB7_SCRIPT}
echo "QC=${QC}" >> ${JOB7_SCRIPT}
echo "ORIG=${FASTQ}" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "DUP=\$QC/duplicate_counts" >> ${JOB7_SCRIPT}
echo "HP=\$QC/homopolymer_files" >> ${JOB7_SCRIPT}
echo "ADAPT=\$QC/qc_adapters/files" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "DUPOUT=\$QC/duplicate_summary.csv" >> ${JOB7_SCRIPT}
echo "HPOUT=\$QC/homopolymer_summary.csv" >> ${JOB7_SCRIPT}
echo "ADAPTOUT=\$QC/adapter_summary.csv" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "### Combine duplicate summaries" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "cd \$DUP" >> ${JOB7_SCRIPT}
echo "FLAG=0" >> ${JOB7_SCRIPT}
echo "for FILE in *.csv" >> ${JOB7_SCRIPT}
echo "do" >> ${JOB7_SCRIPT}
echo "      if [ \$FLAG == 0 ]" >> ${JOB7_SCRIPT}
echo "      then" >> ${JOB7_SCRIPT}
echo '          cat $FILE | sed -e "s/_duplicate_counts_summary.csv//g" > $DUPOUT' >> ${JOB7_SCRIPT}
echo "          FLAG=1" >> ${JOB7_SCRIPT}
echo "      else" >> ${JOB7_SCRIPT}
echo '          tail -n +2 $FILE | sed -e "s/_duplicate_counts_summary.csv//g" >> $DUPOUT' >> ${JOB7_SCRIPT}
echo "      fi" >> ${JOB7_SCRIPT}
echo "done" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "### Combine homopolymer summaries" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo 'ORIGPATH=$( echo \$ORIG/ | sed -e "s/\//\\\\\//g")' >> ${JOB7_SCRIPT}
echo "cd \$HP" >> ${JOB7_SCRIPT}
echo "FLAG=0" >> ${JOB7_SCRIPT}
echo "for FILE in *.csv" >> ${JOB7_SCRIPT}
echo "do" >> ${JOB7_SCRIPT}
echo "      if [ \$FLAG == 0 ]" >> ${JOB7_SCRIPT}
echo "      then" >> ${JOB7_SCRIPT}
echo '          cat $FILE | sed -e "s/${ORIGPATH}//g" > $HPOUT' >> ${JOB7_SCRIPT}
echo "          FLAG=1" >> ${JOB7_SCRIPT}
echo "      else" >> ${JOB7_SCRIPT}
echo '          tail -n +2 $FILE | sed -e "s/${ORIGPATH}//g" >> $HPOUT' >> ${JOB7_SCRIPT}
echo "      fi" >> ${JOB7_SCRIPT}
echo "done" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "### Combine adapter summaries" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo "cd \$ADAPT" >> ${JOB7_SCRIPT}
echo "FLAG=0" >> ${JOB7_SCRIPT}
echo "for FILE in *.csv" >> ${JOB7_SCRIPT}
echo "do" >> ${JOB7_SCRIPT}
echo "      if [ \$FLAG == 0 ]" >> ${JOB7_SCRIPT}
echo "      then" >> ${JOB7_SCRIPT}
echo '          cat $FILE > $ADAPTOUT' >> ${JOB7_SCRIPT}
echo "          FLAG=1" >> ${JOB7_SCRIPT}
echo "      else" >> ${JOB7_SCRIPT}
echo '          tail -n +2 $FILE >> $ADAPTOUT' >> ${JOB7_SCRIPT}
echo "      fi" >> ${JOB7_SCRIPT}
echo "done" >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}
echo >> ${JOB7_SCRIPT}


###### Make QC plots

JOB8_NAME=qc_plots
JOB8_TIME=2:00:00
JOB8_NODES=1
JOB8_MEM=8gb
JOB8_SCRIPT=$SCRIPTS/run_make_qc_plots.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB8_NAME}/g" | 
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB8_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB8_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB8_MEM}/g" > ${JOB8_SCRIPT}
echo "#SBATCH --array=1" >> ${JOB8_SCRIPT}
echo "#" >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo 'module load python/2.7.6' >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo "## Set directories" >> ${JOB8_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB8_SCRIPT}
echo "QC=${QC}" >> ${JOB8_SCRIPT}
echo "PROJNAME=\"${PROJNAME}\"" >> ${JOB8_SCRIPT}
echo "OUTPUT=\$QC/qc_plots" >> ${JOB8_SCRIPT}
echo "if [ ! -e \$OUTPUT ]; then mkdir -p \$OUTPUT; fi" >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo "SCRIPTS=$QCPIPE/scripts" >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo "## Set summary files" >> ${JOB8_SCRIPT}
echo "ADAPT=\$QC/adapter_summary.csv" >> ${JOB8_SCRIPT}
echo "HP=\$QC/homopolymer_summary.csv" >> ${JOB8_SCRIPT}
echo "DUPS=\$QC/duplicate_summary.csv" >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo "## Update sample IDs in adapter and homopolymer summaries" >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo "DESIGN_FILE=\$PROJ/design_files/qc_design_file.csv" >> ${JOB8_SCRIPT}
echo 'while read p; do' >> ${JOB8_SCRIPT}
echo '   IFS="," read -ra ARRAY <<< "$p"' >> ${JOB8_SCRIPT}
echo '   FILE=${ARRAY[1]}' >> ${JOB8_SCRIPT}
echo '   NAME=${ARRAY[2]}' >> ${JOB8_SCRIPT}
echo '   sed -i "s/${FILE}/${NAME}/g" $ADAPT' >> ${JOB8_SCRIPT}
echo '   sed -i "s/${FILE}/${NAME}/g" $HP' >> ${JOB8_SCRIPT}
echo 'done < ${DESIGN_FILE}' >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo '### Make plots for duplicates' >> ${JOB8_SCRIPT}
echo 'python $SCRIPTS/duplicates_plots.py --input $DUPS --obar $OUTPUT/duplicates_by_sample.png --ohist $OUTPUT/duplicates_distribution.png --name "$PROJNAME" ' >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo '### Make plots for adapters ' >> ${JOB8_SCRIPT}
echo 'python $SCRIPTS/adapter_plots.py --input $ADAPT --obar $OUTPUT/adapter_content_by_sample.png --odensity $OUTPUT/adapter_content_distribution.png --name "$PROJNAME" ' >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
echo '### Make plots for homopolymers ' >> ${JOB8_SCRIPT}
echo 'python $SCRIPTS/homopolymer_plots.py --input $HP --obar $OUTPUT/homopolymer_content_by_sample.png --odensity $OUTPUT/homopolymer_content_distribution.png --name "$PROJNAME" ' >> ${JOB8_SCRIPT}
echo >> ${JOB8_SCRIPT}
