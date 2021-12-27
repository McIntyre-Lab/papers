#!/bin/sh
################################################################################################
#
#  Title: MCINTYRE LAB FASTQ QC PIPELINE (NO ERCC CONTROLS), SLURM VERSION
#  Author: Jeremy R. B. Newman
#
#  Description: This script will automatically generate abd submit the SBATCH scripts for running
#  QC on short read Illumina sequencing data. It takes Gzipped FASTQ files and uncompresses them
#  (see below if starting with uncompressed (non-Gzipped) FASTQ files), then analyzes them with
#  FASTQC and custom Python scripts to assess the quality of the data. It also aligns the reads
#  to ERCC spike-in control sequences (single ended alignments) and generates additional plots
#  for the user to assess the overall quality of the data.
#
#  NOTE: This version of the script will NOT align reads to ERCC spike-in control sequences (mix 1). 
#  If you want to run ERCC alignments, you should review and run the run_qc_w_ercc.sh script.
#
#################################################################################################

## Set QC pipeline folder
QCPIPE=/ufrc/mcintyre/share/qc_pipeline

echo "Setting variables"
#################################################################################################
#  SECTION 1: VARIABLES TO BE SET BY THE USER
#
#  In this section, the user will need to specify the value of several variables required for
#  the scripts to run correctly, including the project folder, original data folder, and where
#  to put QC outputs. Prior to running the script, the user will also need to generate a
#  "master" design file for this script to used to generate a QC design file. See SECTION 2
#  more details.
#
#################################################################################################

# Set user email address
EMAIL=nml@ufl.edu

# Set the project folder
PROJ=/ufrc/mcintyre/share/etoh_srna

# Set location of original data. This is the location containing the FASTQ/FASTQ.GZ files
ORIG=$PROJ/original_data/all_fq_files

# Set folder for putting Un-Gzipped FASTQs
# NOTE: if you are starting with FASTQ files that are NOT compressed/gzipped, set this to the
# original data folder ($ORIG)

FASTQ=$PROJ/unzipped_chipSeq_FQ
    if [ ! -e $FASTQ ]; then mkdir -p $FASTQ; fi

# Set folder for QC outputs
QC=$PROJ/qc_chipSeq
    if [ ! -e $QC ]; then mkdir -p $QC; fi

# Set folder for SLURM job logs
SLURMLOGSDIR=$PROJ/scripts/SLURM_LOGS
    if [ ! -e $SLURMLOGSDIR ]; then mkdir -p $SLURMLOGSDIR; fi

# Set the PBS array. As the QC design file has a header row, this should start at 2
# and finished at N+1, where N is the number of samples to run
JOBARRAYIDS=2-73

# Set the project name. Keep this in double quotes (Python limitation)
PROJNAME="etoh_srna_chip"

#################################################################################################
#  SECTION 2: DESIGN FILE GENERATION
#
#  This QC pipeline will use a design file containing two variables: the filename of the FASTQ
#  file (e.g. Sample1_P1_WA01_ATGATG_R1_001.fastq) and a name/label for the file. This name/label
#  can be the filename or something more descriptive (e.g. "P1_WA01_Dmel_ETOH_Rep1_Read1").
#  The user will need to first make a design file that contains the variables they want to use
#  for the filename and name/label, and specify how this script should use this information to
#  build the QC design file. These labels should be unique to each FASTQ file.
#  (Hint: for paired-end data, include the read pair ("R1", "R2") in the label variable)
#
#  An example:
#
#  User-created design file has the following columns:
#      filename, plate, well, sample_id, treatment, replicate, read_pair
#
#  In the subsection "USER-DEFINED DESIGN VARIABLES" (see below) the user will do the following:
#      FILENAME=${ARRAY[0]} 
#      PLATE=${ARRAY[1]}
#      WELL=${ARRAY[2]}
#      SAMPLE=${ARRAY[3]}
#      TREATMENT=${ARRAY[4]}
#      REP=${ARRAY[5]}
#      READPAIR=${ARRAY[6]}


#
#  The user will then define the variables ${FASTQ_FILE} and ${FASTQ_LABEL}, which are used by
#  this script to make the QC design file:
#
#  FASTQ_FILE=${FILENAME}
#  FASTQ_LABEL=${PLATE}_${WELL}_${SAMPLE}_${TREATMENT}_${REP}_${READPAIR}
#
#  The script will then generate the QC design file.
#
################################################################################################
echo "Creating QC design file"

# Location of the user-create design file
MAIN_DESIGN=$PROJ/design_files/DF_trt_and_fq_info.csv

# Run DOS2UNIX on user-created design file to be sure
dos2unix ${MAIN_DESIGN}

## Location of QC design file.
## WARNING: DO NOT CHANGE THIS OR SCRIPTS WILL BREAK.
QC_DESIGN=$PROJ/design_files/qc_design_file.csv

# Header row for the QC design file
echo "filename,label" > ${QC_DESIGN}

# Iterate through the user-created design file and append to QC design file.
# You should change the numbers in {2..97} to match the minimum and maximum numbers
# in the JOBARRAYIDS from Section 1.

for ARRAYID in {74..145}
do
    DESIGN=$(cat ${MAIN_DESIGN} | head -n ${ARRAYID} | tail -n1 )
    IFS=',' read -ra ARRAY <<< "$DESIGN"

# USER-DEFINED DESIGN VARIABLES (you will likely need to change this!)
    FQ=${ARRAY[7]}

    FASTQ_FILE=${FQ}
    FASTQ_LABEL=${FQ}

    echo ${FASTQ_FILE},${FASTQ_LABEL} >> ${QC_DESIGN}
done

#################################################################################################
#  SECTION 3: SCRIPT SUBMISSION
#
#  This section runs the scripts to build your SBATCHs and submit them to the scheduler.
#
#  WARNING: This section should not require any alterations. If a user wishes to edit this
#  section, they do so at their own peril :)
#
#################################################################################################

# Convert PBSLOGSDIR into a SED-friendly format (do not change)
LOGDIR=$( echo $SLURMLOGSDIR | sed -e 's/\//\\\//g')

echo "Building scripts"
### Build scripts
source $QCPIPE/scripts/build_scripts.sh


echo "Submitting SBATCHs"
# Submit generated SBATCHs
### NOTE: the sleep commands here are a fail-safe in case the scheduler takes more than a second to add a job to the queue
### and is there to prevent the next job to read the wrong JOB_ID

#### CHECK: If $FASTQ=$ORIG then skip the gunzip!

if [ $FASTQ != $ORIG ]; then
         # Un-gzip FASTQ files and run FASTQC
         JOB0_ID=$(sbatch ${JOB0_SCRIPT})
         echo "${JOB0_SCRIPT} submitted: Job name is ${JOB0_NAME} and the JobID is ${JOB0_ID}"
         sleep 2
         JOB1_ID=$(sbatch --dependency=afterany:${JOB0_ID/Submitted batch job /} ${JOB1_SCRIPT})
         echo "${JOB1_SCRIPT} submitted: Job name is ${JOB1_NAME} and the JobID is ${JOB1_ID}"
         sleep 2
     else
         # Only run FASTQC
         JOB1_ID=$(sbatch ${JOB1_SCRIPT})
         echo "${JOB1_SCRIPT} submitted: Job name is ${JOB1_NAME} and the JobID is ${JOB1_ID}"
         rm ${JOB0_SCRIPT}
         sleep 2
fi

# FASTQC Summaries
JOB2_ID=$(sbatch --dependency=afterany:${JOB1_ID/Submitted batch job /} $JOB2_SCRIPT)
echo "${JOB2_SCRIPT} submitted: Job name is ${JOB2_NAME} and the JobID is ${JOB2_ID}"
sleep 2

# Count duplicates
JOB3_ID=$(sbatch --dependency=afterany:${JOB2_ID/Submitted batch job /} ${JOB3_SCRIPT})
echo "${JOB3_SCRIPT} submitted: Job name is ${JOB3_NAME} and the JobID is ${JOB3_ID}"
sleep 2

# Identify homopolymers
JOB4_ID=$(sbatch --dependency=afterany:${JOB3_ID/Submitted batch job /} ${JOB4_SCRIPT})
echo "${JOB4_SCRIPT} submitted: Job name is ${JOB4_NAME} and the JobID is ${JOB4_ID}"
sleep 2

# FASTQ to FASTA
JOB5_ID=$(sbatch --dependency=afterany:${JOB4_ID/Submitted batch job /} ${JOB5_SCRIPT})
echo "${JOB5_SCRIPT} submitted: Job name is ${JOB5_NAME} and the JobID is ${JOB5_ID}"
sleep 2

# BLAT Illumina adapters
JOB6_ID=$(sbatch --dependency=afterany:${JOB5_ID/Submitted batch job /} ${JOB6_SCRIPT})
echo "${JOB6_SCRIPT} submitted: Job name is ${JOB6_NAME} and the JobID is ${JOB6_ID}"
sleep 2

# Summarize QC files
JOB7_ID=$(sbatch --dependency=afterany:${JOB6_ID/Submitted batch job /} ${JOB7_SCRIPT})
echo "${JOB7_SCRIPT} submitted: Job name is ${JOB7_NAME} and the JobID is ${JOB7_ID}"
sleep 2

# Make QC plots
JOB8_ID=$(sbatch --dependency=afterany:${JOB7_ID/Submitted batch job /} ${JOB8_SCRIPT})
echo "${JOB8_SCRIPT} submitted: Job name is ${JOB8_NAME} and the JobID is ${JOB8_ID}"
sleep 2

# Delete ERCC qsubs generated
rm ${JOB9_SCRIPT}
rm ${JOB10_SCRIPT}
rm ${JOB11_SCRIPT}
rm ${JOB12_SCRIPT}

### CLEAN UP
echo "Clean-up script has been generated but not submitted. It is located at ${JOB13_SCRIPT}"
echo "The user should verify that all previous QC jobs were successfully run before submitting"
echo "the clean-up script to the scheduler."
sleep 3

echo "All QC jobs generated/submitted!"
