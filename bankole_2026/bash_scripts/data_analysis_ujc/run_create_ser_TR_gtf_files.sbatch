#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=koppGTF
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/koppGTF_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=3-49
##SBATCH --array=2

## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
ALN=${PROJ}/aln_ser_dorReads_2_dser1

## Get design file information (2 - 49)
DESIGN_FILE=$PROJ/design_files/design_file_sample2BC_w_origDataPath_04amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

SAM=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.sam
BAM=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.bam
GTF=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.gtf
BED=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.bed

module purge
module load samtools/1.10 bedtools/2.29.2

## SAM to BAM
samtools view -b -F 2048 ${SAM} > ${BAM}

## BAM to BED
bedtools bamtobed -split -i ${BAM} > ${BED}

## BED to GTF
SOURCE=${SAMPLEID}
awk -v source=${SOURCE} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

## Remove intermediate files
rm ${BED} ${BAM}
