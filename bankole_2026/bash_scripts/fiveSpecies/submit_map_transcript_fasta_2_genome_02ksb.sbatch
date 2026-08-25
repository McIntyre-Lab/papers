#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=map_2_one
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/map_2_one_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb
#SBATCH --account=mcintyre
#SBATCH --array=2-6

## Description: Map transcript FASTA to genome FASTA. SAM -> BAM -> BED -> GTF. Bin unmapped transcripts.

module purge
ml samtools/1.18 bedtools/2.30.0
module load minimap2/2.24

#export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
SCRIPTS=$PROJ/scripts
REFERENCES=/blue/mcintyre/share/references

# loop through all species and self-map -> SAM/BAM/BED/GTF -> id_ujc
# count: number mapping, number UJC

DESIGN_FILE=$PROJ/design_files/df_map_transcripts_to_genome_03ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

OUTNAME=${ARRAY[0]}
GENOME_PATH=${ARRAY[1]}
TR_FASTA=${ARRAY[2]}

# References Folder for genome
OUTD=${ARRAY[3]}

ROZ=$OUTD/roz_mapping_${OUTNAME}
	mkdir -p $ROZ

ANNO_PATH=${ARRAY[4]}


echo "Aligning $TR_FASTA to $GENOME_PATH"
# -a output sam format, -x splice for longread spliced alignment, removed --for-only that only maps to forward strand of ref, no secondary alignments, -C 5 cost for non-canonical splicing
minimap2 \
    -t 12 \
    -a \
    -x splice \
    --secondary=no \
    -C 5 \
    ${GENOME_PATH} \
    ${TR_FASTA} \
    > ${ROZ}/${OUTNAME}_unfilt_minimap2.sam 2>$ROZ/${OUTNAME}_unfilt_minimap2.log

echo "Binning unmapped"
# Bin unmapped transcripts
SAM=$ROZ/${OUTNAME}_unfilt_minimap2.sam
samtools view -f 4 $SAM | awk '{print $1}' | sort | uniq > $OUTD/list_${OUTNAME}_unmapped.txt

echo "Num unmapped transcripts:"
wc -l $OUTD/list_${OUTNAME}_unmapped.txt

#grep -f $ROZ/${OUTNAME}_unmapped_transcripts.txt $ANNO_PATH > $OUTD/${OUTNAME}_unmapped.gtf

echo "Converting SAM to GTF"
# SAM -> BAM -> BED -> GTF
samtools view -b -F 2048 $SAM > ${ROZ}/${OUTNAME}_minimap2.bam
bedtools bamtobed -split -i ${ROZ}/${OUTNAME}_minimap2.bam > ${ROZ}/${OUTNAME}_minimap2.bed
BED=${ROZ}/${OUTNAME}_minimap2.bed
GTF=${OUTD}/${OUTNAME}_roz.gtf

#chrom,source,exon,start,end,.,strand,.,geneID;transcriptID
awk -v source=${OUTNAME} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

echo "Num mapped transcripts:"
awk '$3~"exon"' $GTF | awk '{print $12}' | sort | uniq | wc -l

rm -r $ROZ
