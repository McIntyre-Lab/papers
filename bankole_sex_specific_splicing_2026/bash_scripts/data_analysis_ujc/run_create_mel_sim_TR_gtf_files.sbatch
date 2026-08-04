#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=rmgGTF
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/rmgGTF_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=3-109
##SBATCH --array=2

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_bc_design_w_origDataPath_04amm.csv
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

if [[ $SAMPLEID == mel* ]]; then
	ALN=${PROJ}/aln_mel_dorReads_2_dmel6
	GENOME=dmel6
else
	ALN=${PROJ}/aln_sim_dorReads_2_dsim2
	GENOME=dsim2
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
