#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=map_2_one_03
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/map_2_one_03_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

## Description: Map ref transcript FASTA to ref genome FASTA. Convert SAM -> BAM -> BED -> GTF. Bin "removed" transcripts.

module purge
ml samtools/1.18 bedtools/2.30.0
module load minimap2/2.24

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
SCRIPTS=$PROJ/scripts
REFERENCES=/blue/mcintyre/share/references

DESIGN_FILE=$PROJ/design_files/df_map_transcripts_to_genome_03ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

OUTNAME=${ARRAY[0]}
GENOME_PATH=${ARRAY[1]}
TR_FASTA=${ARRAY[2]}

# OUTD=References Folder for genome
OUTD=${ARRAY[3]}

ROZ=$OUTD/roz_mapping_${OUTNAME}
	mkdir -p $ROZ

ANNO_PATH=${ARRAY[4]}

echo "Aligning $TR_FASTA to $GENOME_PATH"
#nadja's minimap settings
minimap2 \
	-t $SLURM_CPUS_ON_NODE \
	-a \
	-x splice \
	--secondary=yes \
	-N 200 \
	-p 0.9 \
	-C 5 \
	${GENOME_PATH} \
	${TR_FASTA} \
	> ${ROZ}/${OUTNAME}_unfilt_minimap2.sam 2>$ROZ/${OUTNAME}_unfilt_minimap2.log

# SAM -> BAM -> BED -> GTF
echo "Converting SAM to GTF"
SAM=${ROZ}/${OUTNAME}_unfilt_minimap2.sam

## remove transcripts with mapping quality of 0 (remove header)
samtools view \
	-h \
	-q 1 \
	${SAM} \
	> ${ROZ}/${OUTNAME}_noMQ0.sam

## sam2bam, exclude reads using flag 2052 (omit unmapped and supplementary)
samtools view \
	-b \
	-F 2052 \
	${ROZ}/${OUTNAME}_noMQ0.sam \
	> ${ROZ}/${OUTNAME}_minimap2.bam

# bam2bed
bedtools bamtobed -split -i ${ROZ}/${OUTNAME}_minimap2.bam > ${ROZ}/${OUTNAME}_minimap2.bed

#bed2gtf
BED=${ROZ}/${OUTNAME}_minimap2.bed
GTF=${OUTD}/${OUTNAME}_noGeneID.gtf

#chrom,source,exon,start,end,.,strand,.,geneID;transcriptID
awk -v source=${OUTNAME} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

# Print Counts
echo "Num mapped transcripts:"
awk '$3~"exon"' $GTF | awk '{print $12}' | sort | uniq | wc -l

echo "Binning removed (unmapped, supplementary, mq=0)"
# Bin removed transcripts
awk '!/^@/ {print $1}' ${ROZ}/${OUTNAME}_unfilt_minimap2.sam | sort | uniq > ${ROZ}/list_full_tr.txt
awk '$3~"exon"' $GTF | awk '{print $12}' | tr -d '";' | sort | uniq > ${ROZ}/list_subset_tr.txt

# Find the transcripts from the SAM file that are not in the GTF and list them
comm -23 ${ROZ}/list_full_tr.txt ${ROZ}/list_subset_tr.txt > $OUTD/list_${OUTNAME}_removed_tr.txt

echo "Num removed transcripts:"
wc -l $OUTD/list_${OUTNAME}_removed_tr.txt

rm -r $ROZ
