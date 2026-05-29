#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=isoRefine
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/isoRefine_%A-%a.out
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-61


### isoseq3 analysis T1D
## refine
## trim polyA tails
## concatemer ID and removal

module load isoseq3/3.8.2
module load samtools/1.12

### Set Project Directory
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio


## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_bc_2_run.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]} 
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}
SAMPLEID=${ARRAY[3]}
BC=${ARRAY[4]}

FULLID=${NAME}_demuxed
SAMPLEID_NOHYPHEN=${SAMPLEID//-/_}

### Set Directories
OUTD=$PROJ/isoseq3_refine_by_sample/${RUN}
    mkdir -p $OUTD
ROZ=${OUTD}/roz
    mkdir -p ${ROZ}
IND=$PROJ/isoseq3_library_processing
IN=$PROJ/original_data


date


if [[ ${RUN} == "first_run" ]] ; then
    echo "run is $RUN"
    PRIMERS=${IN}/${RUN}/${RUNID}/IsoSeq_Primers_12_Barcodes_v1.fasta
else
    echo "run is $RUN"
    PRIMERS=${IN}/${RUN}/${RUNID}/barcodes.fa
fi


##Combine 5p and 3p bam files

if [[ ${RUN} == "first_run" ]] ; then
     samtools merge ${ROZ}/${FULLID}.${SAMPLEID}.bam $IND/${RUN}/demultiplexed/${FULLID}.${BC}*.bam
else
     samtools merge ${ROZ}/${FULLID}.${SAMPLEID}.bam $IND/${RUN}/demultiplexed/${FULLID}.${SAMPLEID_NOHYPHEN}*.bam
fi




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
    ${ROZ}/${FULLID}.${SAMPLEID}.bam \
    $PRIMERS \
    $OUTD/${SAMPLEID}_${RUN}.flnc.bam

#rm ${ROZ}/${FULLID}.${SAMPLEID}.bam 
