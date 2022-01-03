#!/bin/bash

###### Build clean-up qsub

JOB13_NAME=qc_cleanup
JOB13_TIME=2:00:00
JOB13_NODES=1
JOB13_MEM=4gb
JOB13_SCRIPT=$SCRIPTS/clean_up_qc.sbatch

less $TEMPLATES/slurm_header.txt |
    sed -e "s/EMAIL_HERE/${EMAIL}/g" |
    sed -e "s/JOBNAME_HERE/${JOB13_NAME}/g" |
    sed -e "s/PBS_LOGS_HERE/${LOGDIR}/g" |
    sed -e "s/WALLTIME_HERE/${JOB13_TIME}/g" |
    sed -e "s/NODES_HERE/${JOB13_NODES}/g" |
    sed -e "s/MEM_HERE/${JOB13_MEM}/g" > ${JOB13_SCRIPT}
echo "#PBS -t 1" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "##### Script to clean up QC folder by removing unneeded intermediate files" >> ${JOB13_SCRIPT}
echo "##### The following files will be kept:" >> ${JOB13_SCRIPT}
echo "## FASTQC: PDFs, pass-fail summary CSV, list of Kmers and overrepresented sequences per sample" >> ${JOB13_SCRIPT}
echo "## QC: adapter, homopolymer and duplicates summaries and plots" >> ${JOB13_SCRIPT}
echo "## ERCC: concentration plots, merged coverage CSV and SAS datasets" >> ${JOB13_SCRIPT}
echo "## FASTQ Split Dups PE (this is in the project folder)" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "##### The following files will be deleted:" >> ${JOB13_SCRIPT}
echo "## Gunzipped FASTQ files" >> ${JOB13_SCRIPT}
echo "## FASTQC: HTML files" >> ${JOB13_SCRIPT}
echo "## QC: All by-sample files, FASTA files, other intermediate files" >> ${JOB13_SCRIPT}
echo "## ERCC: alignments, mpileups, coverage counts per sample" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "## Set directories" >> ${JOB13_SCRIPT}
echo "PROJ=${PROJ}" >> ${JOB13_SCRIPT}
echo "QC=${QC}" >> ${JOB13_SCRIPT}
echo "FASTQ=${FASTQ}" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "## Remove GUnzipped FASTQ files" >> ${JOB13_SCRIPT}
echo "rm -r \${FASTQ}" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "## Remove FASTQC intermediate files" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/fastqc/fastqc_files" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "## Remove QC intermediate files" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "# Duplicate counts" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/fastq_split_dups" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/duplicate_counts" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "# Homopolymers" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/homopolymer_files" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/homopolymer_logs" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "# Adapters" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/fasta_reads" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/qc_adapters" >> ${JOB13_SCRIPT}
echo >> ${JOB13_SCRIPT}
echo "## Remove ERCC intermediate files" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/ercc/aln_ERCC_all_se" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/ercc/ercc_se_mpileups" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/ercc/bam_files_ercc" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/ercc/coverage_counts_ercc" >> ${JOB13_SCRIPT}
echo "rm -r \$QC/sas_temp" >> ${JOB13_SCRIPT}

