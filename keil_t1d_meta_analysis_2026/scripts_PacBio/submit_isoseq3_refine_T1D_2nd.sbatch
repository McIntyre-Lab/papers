#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=isoRefine
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/isoRefine_%A-%a.out
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-13


### isoseq3 analysis T1D -- Second Run!!
## refine
## trim polyA tails
## concatemer ID and removal

module load isoseq3/3.4.0

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
OUTR=$PROJ/isoseq3_refine/second_run
    mkdir -p $OUTR

IND=$PROJ/isoseq3_library_processing/second_run/demultiplexed

PRIMERS=$PROJ/original_data/second_run/4470/barcodes.fa

## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_barcode_2nd_3rd.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
BC=${ARRAY[2]}

date

echo "
*** Refine - remove polyA and concatemers***
"
## Refine
## I believe samples have polyA so use --require-polya parameter
## using default parameters:
        ##  --min-polya-length  Minimum poly(A) tail length. [20]
        ##  --min-rq            Minimum CCS RQ. Default is -1, deactivated [-1]

isoseq3 refine \
    --num-threads 24 \
    --require-polya \
    $IND/m64047_210212_140334_demuxed.${BC}.bam	\
    $PRIMERS \
    $OUTR/${SAMPLEID}_m64047_210212_140334_demuxed.${BC}.flnc.bam

