#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=tofu2
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/tofu_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=4
#SBATCH --array=1-11

### Tofu2 Cluster and Alignment Processing
### all samples vs mo17 cau reference

module load samtools tofu2/20180629 

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/tofu2_mo17_cau
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IND=$PROJ/mapping_minimap2_mo17_cau

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

REF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa

## Set input and output directories to the correct species
OUT=$OUTD/${SAMPLE_NAME}
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IN=$IND/$GENO/$TRT

date
echo "***Combine Polish BAM***
"
## Combine Polish bam output files
## all read names are staying in the form "transcript/#"

if [ ! -e $OUT/${SAMPLE_NAME}.polished.all.bam ]; then
    samtools merge $OUT/${SAMPLE_NAME}.polished.all.bam \
        $(ls $PROJ/isoseq3_analysis/polish_by_individual/$GENO/$TRT/${SAMPLE_NAME}.polished.transcriptset.?.bam \
        $PROJ/isoseq3_analysis/polish_by_individual/$GENO/$TRT/${SAMPLE_NAME}.polished.transcriptset.??.bam)
    echo "    Creating combined polish bam file"
else
    echo "    Combined polish bam file already exists"
fi

date
echo "***Make Cluster Report***
"
## Make Cluster Report of Polished Reads
## input combined polished file to make cluster_report.csv
## tofu2 script issue with opening bam file - requires "check_sq=False" in pysam.Samfile(polished_bam)
## added an output argument to isoseq3_cluster_report_cup.py

python $PROJ/scripts/pacbio/isoseq3_make_cluster_report_cup.py $OUT/${SAMPLE_NAME}.polished.all.bam $OUT/${SAMPLE_NAME}

date
echo "***Sort Mapped SAM***
"
## Sort mapped SAM for further processing
## Remove header lines with samtools view and remove all secondary alignments with -F 2048
## Unmapped reads are still included in the sorted sam output file

samtools view -F 2048 $IN/${SAMPLE_NAME}.polished.all.hq.mapped.sam | sort -k3,3 -k4,4n > $IN/${SAMPLE_NAME}.polished.all.hq.mapped.sorted.sam

date
echo "***Collapse Isoforms***
"
## Collapse Isoforms From Sorted Sam
## requires the (unzipped) fastq from before alignment and the mapped sorted sam from after alignment
## --dun-merge-5-shorter collapses shorter 5' transcripts (selects for longest)
## outputs a gff file that is actually in gtf format

collapse_isoforms_by_sam.py \
    --input $IN/${SAMPLE_NAME}.polished.all.hq.fastq \
    --fq \
    -s $IN/${SAMPLE_NAME}.polished.all.hq.mapped.sorted.sam \
    --dun-merge-5-shorter \
    -o $OUT/${SAMPLE_NAME}

date
echo "***Get Abundance***
"
## Get Abundance of Collapsed Reads
## Use cluster report to get abundance information of collapsed reads

get_abundance_post_collapse.py $OUT/${SAMPLE_NAME}.collapsed $OUT/${SAMPLE_NAME}.cluster_report.csv

date
