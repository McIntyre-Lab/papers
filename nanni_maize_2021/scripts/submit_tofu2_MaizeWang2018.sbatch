#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=tofu2
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/tofu2_MaizeWang2018_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=4
#SBATCH --array=1-5

### Tofu2 Cluster and Alignment Processing
### 9 samples total

module purge
module load samtools/1.10 tofu2/20180629

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval
POLISHOUT=${OUTD}/isoseq3_library_processing/polish
MAPOUT=${OUTD}/minimap2
TOFUOUT=${OUTD}/tofu2
    mkdir -p ${TOFUOUT}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Get reference genome FASTA
if [[ ${SPECIES} == "maize" ]]; then
    REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
else
    REF=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
fi

date

## Combine Polish bam output files
## all read names are staying in the form "transcript/#"

if [ ! -e ${POLISHOUT}/${SAMPLEID}.polished.all.bam ]; then
    samtools merge ${POLISHOUT}/${SAMPLEID}.polished.all.bam \
        $(ls ${POLISHOUT}/${SAMPLEID}.polished.transcriptset.?.bam \
        ${POLISHOUT}/${SAMPLEID}.polished.transcriptset.??.bam)
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

python $PROJ/scripts/pacbio/isoseq3_make_cluster_report_cup.py \
    ${POLISHOUT}/${SAMPLEID}.polished.all.bam \
    ${TOFUOUT}/${SAMPLEID}
date

echo "***Sort Mapped SAM***
"
## Sort mapped SAM for further processing
## Remove header lines with samtools view and remove all secondary alignments with -F 2048
## Unmapped reads are still included in the sorted sam output file

samtools view -F 2048 ${MAPOUT}/${SAMPLEID}.polished.all.hq.mapped.sam | \
    sort -k3,3 -k4,4n \
    > ${MAPOUT}/${SAMPLEID}.polished.all.hq.mapped.sorted.sam
date

echo "***Collapse Isoforms***
"
## Collapse Isoforms From Sorted Sam
## requires the (unzipped) fastq from before alignment and the mapped sorted sam from after alignment
## --dun-merge-5-shorter collapses shorter 5' transcripts (selects for longest)
## outputs a gff file that is actually in gtf format

collapse_isoforms_by_sam.py \
    --input ${MAPOUT}/${SAMPLEID}.polished.all.hq.fastq \
    --fq \
    -s ${MAPOUT}/${SAMPLEID}.polished.all.hq.mapped.sorted.sam \
    --dun-merge-5-shorter \
    -o ${TOFUOUT}/${SAMPLEID}

date
echo "***Get Abundance***
"
## Get Abundance of Collapsed Reads
## Use cluster report to get abundance information of collapsed reads

get_abundance_post_collapse.py \
    ${TOFUOUT}/${SAMPLEID}.collapsed \
    ${TOFUOUT}/${SAMPLEID}.cluster_report.csv
date
