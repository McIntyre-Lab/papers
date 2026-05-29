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

### isoseq3 analysis - library processing
### have 1 subreads.bam file that should be barcoded. 3 cell types, unknown # of samples
### starting input is a *.subreads.bam file - sequences include subreads separated by primer and adapter sequences
### lima will demultiplex


module load isoseq3/3.4.0 pbccs/5.0.0 lima/2.0.0 

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
OUTD=$PROJ/isoseq3_library_processing
    mkdir -p $OUTD
OUTL=$OUTD/demultiplexed
    mkdir -p $OUTL
IN=$PROJ/original_data/3988
ROZ=$PROJ/roz_isoseq3
    mkdir -p $ROZ

PRIMERS=$IN/IsoSeq_Primers_12_Barcodes_v1.fasta

LOG=$OUTD/log
    mkdir -p $LOG

#date

echo"
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
    --log-file $LOG/m64047_200918_062828.log \
    --report-file $OUTD/m64047_200918_062828_report.txt \
    $IN/m64047_200918_062828.subreads.bam \
    $OUTD/m64047_200918_062828.subreads.ccs.bam


date

echo "
***Primer Removal and Demultiplexing***

Input: $OUTD/m64047_200918_062828.subreads.ccs.bam
"
## --isoseq activate specialized isoseq mode for primer removal and demultiplexing
## --different Only keep different barcodes in a pair in output. Enforces --min-passes ≥ 1
## --split-names split output by resolved barcode pair name
## min-passes 1

lima --num-threads 12 \
    --isoseq \
    --different \
    --peek-guess \
    $OUTD/m64047_200918_062828.subreads.ccs.bam \
    $PRIMERS \
    $OUTL/m64047_200918_062828_demuxed.bam

rm -rf $ROZ
