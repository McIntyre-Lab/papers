#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=bcTO
#SBATCH --account=mcintyre-cid
#SBATCH --qos=mcintyre-cid
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/bcTO_%A_%a.out
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24gb
#SBATCH --partition=gpu
#SBATCH --gpus=a100:2
#SBATCH --array=3
date;hostname;pwd

#rebasecalling

module load pod5/0.1.5 cuda pytorch dorado/0.5.2 guppy/6.4.2

PROJ=/blue/mcintyre/share/transcript_ortholog
OUT=$PROJ/dorado_basecalling

## download models - uses pod5 files to automatically select models
    ## analyte type: dna
    ## pore type: for min106D r9.4.1
    ## chemistry type: kit used kit 9 chemistry so need e8 chemistry type
    ## model type: fast = quickest, sup = most accurate, hac = balance between speed and accuracy
    ## dna_r9.4.1_e8_hac@v3.6
#dorado download --directory $OUT/dorado_models/ --data $POD5//


## design rows 
DESIGN_FILE=$PROJ/design_files/list_ont_dirs_flagPod5.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
TR=${ARRAY[1]}
FC=${ARRAY[2]}
DATE=${ARRAY[3]}
DIR=${ARRAY[4]}

date
echo "starting basecalling on pod5 files from ${FC}_${DATE}_TR${TR}
"
## model options:
    # dna_r9.4.1_e8_fast@v3.4 = quickest,
    # dna_r9.4.1_e8_sup@v3.6 = most accurate,
    # dna_r9.4.1_e8_hac@v3.3  = balance between speed and accuracy
MODEL=/blue/mcintyre/share/pxs_lmm_dros_data/dorado_basecalling/dorado_models/dna_r9.4.1_e8_sup@v3.6

    ## model and data dir are positional args
    ## --min-qscore  Discard reads with mean Q-score below this threshold. [default: 0]
    ## --min-qscore  default 0
    ## --kit-name    enable barcoding with kit
    ## --no-trim     Skip trimming of barcodes, adapters, and primers.
    ## --resume-from if crashes, rename incomplete bam and restart

mkdir -p $OUT/basecalls_1to1

    ## basecall on orig pod5 files
date
   if [[ $FLAG == 0 ]]; then
    dorado basecaller \
        --recursive \
        --device "cuda:0,1" \
        --kit-name SQK-PCB109 \
        --trim none \
        ${MODEL} \
        $OUT/pod5_files_1to1/$DIR \
        > $OUT/basecalls_1to1/${FC}_${DATE}_TR${TR}_noTrim.bam
    else
        dorado basecaller \
        --recursive \
        --device "cuda:0,1" \
        --kit-name SQK-PCB109 \
        --trim none \
        ${MODEL} \
        $PROJ/original_data/$DIR \
        > $OUT/basecalls_1to1/${FC}_${DATE}_TR${TR}_noTrim.bam
    fi

date
echo "summarizing basecalls
"
    dorado summary \
        $OUT/basecalls_1to1/${FC}_${DATE}_TR${TR}_noTrim.bam \
        > $OUT/basecalls_1to1/summary_noTrim_${FC}_${DATE}_TR${TR}.txt


date
echo "demuxing ${FC}_${DATE}_TR${TR}
"
    # --no-classify skipe barcode classification, demux based on existing classification in reads
    # --kit-name SQK-PCB109 : cannot be used with --no-classify
    # --no-trim : default is barcode trimming

## testing adding --kit-name and removing --no-classify to demux
mkdir -p $OUT/demux_1to1/${FC}_${DATE}_TR${TR}

    dorado demux \
        $OUT/basecalls_1to1/${FC}_${DATE}_TR${TR}_noTrim.bam \
        --output-dir $OUT/demux_1to1/${FC}_${DATE}_TR${TR} \
        --no-classify \
        --emit-fastq

## count in summary
#filename
#read_id
#run_id
#channel
#mux
#start_time
#duration
#template_start
#template_duration
#sequence_length_template
#mean_qscore_template
#barcode
