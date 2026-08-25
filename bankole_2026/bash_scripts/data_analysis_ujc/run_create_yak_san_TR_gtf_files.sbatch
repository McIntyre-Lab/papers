#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=rlrGTF
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/rlrGTF_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
##SBATCH --array=2-9
#SBATCH --array=2

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog

ml samtools/1.10  bedtools/2.29.2

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
SAMPLE=${ARRAY[1]}
BC=${ARRAY[2]}
BCNUM=${ARRAY[3]}
TR=${ARRAY[4]}
TECH=${ARRAY[5]}
DATE=${ARRAY[6]}
ONT_DIR=${ARRAY[7]}

SAMPLEID=${SAMPLE}_TR${TR}

ALN=${PROJ}/aln_dorReads_2_dyak2_or_dsan1

if [[ $SAMPLEID == dsan* ]] ; then
	GENOME=dsan1
else
	GENOME=dyak2
fi

SAM=${ALN}/${SAMPLEID}_2_${GENOME}_noMQ0_uniq.sam
BAM=${ALN}/${SAMPLEID}_2_${GENOME}_noMQ0_uniq.bam
GTF=${ALN}/${SAMPLEID}_2_${GENOME}_noMQ0_uniq.gtf
BED=${ALN}/${SAMPLEID}_2_${GENOME}_noMQ0_uniq.bed

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
