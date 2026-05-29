#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=lima
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/lima_%A.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --cpus-per-task=12
#SBATCH --array=2

### isoseq3 analysis - library processing
### for each forth and fifth runs:

    ### have 1 subreads.bam file that should be barcoded. 3 cell types, unknown # of samples
    ### starting input is a *.subreads.bam file - sequences include subreads separated by primer and adapter sequences
    ### lima will demultiplex


module load isoseq3/3.4.0 pbccs/5.0.0 lima/2.0.0 

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IN=$PROJ/original_data

## design file for 4th and 5th runs - has header!!
DESIGN_FILE=$PROJ/design_files/design_PB_4th_5th_runs.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]}
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}


OUTD=$PROJ/isoseq3_library_processing/${RUN}
    mkdir -p $OUTD
OUTL=$OUTD/demultiplexed
    mkdir -p $OUTL
ROZ=$PROJ/roz_isoseq3/${RUN}
    mkdir -p $ROZ

SAMPLE_PATH=${RUN}/${RUNID}
PRIMERS=$IN/${SAMPLE_PATH}/barcodes.fa

LOG=$OUTD/log
    mkdir -p $LOG

#date

echo "
***Consensus Calling***

Input: ${SUBREADS_BAM}
"
## Consensus Calling
## outputs consensus sequence from each ZMW with at least one full pass in ccs.bam
## --skip-polish Only output the initial draft templates
## --minPasses minimum number of full passes per read (each primer is seen at least once)
## default values:
	## --min-snr 2.5 
	## --min-rq 0.99 (AMM thinks this will be Q20)
	## --min-length 10
	## --max-length 50000

ccs --min-passes 1 \
    --num-threads 12 \
    --log-file $LOG/${NAME}.log \
    --report-file $OUTD/${NAME}_report.txt \
    $IN/${SAMPLE_PATH}/${NAME}.subreads.bam \
    $OUTD/${NAME}.subreads.ccs.bam


date

echo "
***Primer Removal and Demultiplexing***

Input: $OUTD/${NAME}.subreads.ccs.bam
"
## --isoseq activate specialized isoseq mode for primer removal and demultiplexing
## --different Only keep different barcodes in a pair in output. Enforces --min-passes ≥ 1
## --split-names split output by resolved barcode pair name
## min-passes 1

lima --num-threads 12 \
    --isoseq \
    --different \
    --peek-guess \
    $OUTD/${NAME}.subreads.ccs.bam \
    $PRIMERS \
    $OUTL/${NAME}_demuxed.bam

rm -rf $ROZ
