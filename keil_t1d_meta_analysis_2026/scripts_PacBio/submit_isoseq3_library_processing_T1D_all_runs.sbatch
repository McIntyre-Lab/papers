#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=lima
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/lima_%A.%a.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --cpus-per-task=12
#SBATCH --array=3-6

### isoseq3 analysis - library processing
### for each second and third runs:

    ### have 1 subreads.bam file that should be barcoded. 3 cell types, unknown # of samples
    ### starting input is a *.subreads.bam file - sequences include subreads separated by primer and adapter sequences
    ### lima will demultiplex


module load isoseq3/3.4.0 pbccs/6.4.0 lima/2.7.1

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IN=$PROJ/original_data

## design file - has header!!
DESIGN_FILE=$PROJ/design_files/design_PB_all_runs.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]}
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}

OUT=$PROJ/isoseq3_library_processing/${RUN}
    mkdir -p $OUT
OUTD=$OUT/demultiplexed
    mkdir -p $OUTD
LOG=$OUT/log
    mkdir -p $LOG

#date

#echo"
#***Consensus Calling***
#
#Input: $IN/${RUN}/${RUNID}/${NAME}.subreads.bam
#"
## Consensus Calling
## outputs consensus sequence from each ZMW with at least one full pass in ccs.bam
## --skip-polish Only output the initial draft templates
## --minPasses minimum number of full passes per read (each primer is seen at least once)
## default values:
	## --min-snr 2.5 
	## --min-rq 0.99 (AMM thinks this will be Q20)
	## --min-length 10
	## --max-length 50000

#ccs --min-passes 1 \
#    --num-threads 12 \
#    --log-file $LOG/${NAME}.log \
#    --report-file $OUT/${NAME}_report.txt \
#    $IN/${RUN}/${RUNID}/${NAME}.subreads.bam \
#    $OUT/${NAME}.subreads.ccs.bam
#echo"
#*** Consensus Calling Done ***
#"

date

if [[ ${RUN} ==	"first_run" ]] ; then 
    echo "run is $RUN"
    PRIMER=${IN}/${RUN}/${RUNID}/IsoSeq_Primers_12_Barcodes_v1.fasta
else
    PRIMER=${IN}/${RUN}/${RUNID}/barcodes.fa
fi

echo "
***Primer Removal and Demultiplexing***

BAM: $OUT/${NAME}.subreads.ccs.bam
Primer: ${PRIMER}

"

## --isoseq activate specialized isoseq mode for primer removal and demultiplexing
## --different Only keep different barcodes in a pair in output. Enforces --min-passes ≥ 1
## --split-names split output by resolved barcode pair name
## min-passes 1
## hifi-preset: --ccs --min-score 80 --min-end-score 50 --min-ref-span 0.75 plus --same and --ignore-missing-adapters for SYMMETRIC-ADAPTERS
    ## amm replaced isoseq mode with hifi-preset

lima --num-threads 12 \
    --split-named \
    --store-unbarcoded \
    --hifi-preset SYMMETRIC-ADAPTERS \
    --peek-guess \
    $OUT/${NAME}.subreads.ccs.bam \
    ${PRIMER} \
    $OUTD/${NAME}_demuxed.bam
