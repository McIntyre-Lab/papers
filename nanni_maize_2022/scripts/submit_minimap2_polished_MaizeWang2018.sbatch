#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=minimap2
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/minimap2_polish_MaizeWang2018_%A_%a.out
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-5

## Map to maize B73 genome with minimap2

module load minimap/2.12

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/minimap2
    mkdir -p ${OUTD}
LOG=${OUTD}/log
    mkdir -p ${LOG}
IND=$PROJ/MaizeWang2018_transcriptome_eval/isoseq3_library_processing/polish

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

date

## Get reference genome FASTA
if [[ ${SPECIES} == "maize" ]]; then
    REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
else
    REF=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
fi

date
echo "***Combine Polish Fasta.gz***
"
## Combine Polish hq.fasta.gz output files
## hq files are the high quality clusted consensus reads with predicted accuracy >= 0.99
if [ ! -e ${OUTD}/${SAMPLEID}.polished.all.hq.fasta ]; then
    zcat ${IND}/${SAMPLEID}.*.?.hq.fasta.gz ${IND}/${SAMPLEID}.*.??.hq.fasta.gz > ${OUTD}/${SAMPLEID}.polished.all.hq.fasta
    echo "    Creating combined polish fasta file"
else
    echo "    Combined polish fasta file already exists"
fi

date
echo "***Combine Polish Fastq.gz***
"
## Combine Polish hq.fastq.gz output files
## Not used in mapping but used in downstream analysis
## hq files are the high quality clusted consensus reads with predicted accuracy
##  >= 0.99

if [ ! -e ${OUTD}/${SAMPLEID}.polished.all.hq.fastq ]; then
    zcat ${IND}/${SAMPLEID}.*.?.hq.fastq.gz ${IND}/${SAMPLEID}.*.??.hq.fastq.gz > ${OUTD}/${SAMPLEID}.polished.all.hq.fastq
    echo "    Creating combined polish fastq file"
else
    echo "    Combined polish fastq file already exists"
fi

date
echo "***Minimap2 Map***
"
## Minimap2 Map
## -a output in SAM format
## -x splice uses long read spliced alignment settings
##     meaning long deletions are considered introns, no long insertions,
##     and gap costs are changed accordingly
## -u f option requires mapping of the forward strand
## --secondary=no does not include secondary alignments
## -C 5 modify cost of non-canonical splice junctions (recommended for Isoseq processed PacBio)

minimap2 \
    -t 12 \
    -a -x splice -u f \
    --secondary=no -C 5 \
    ${REF} ${OUTD}/${SAMPLEID}.polished.all.hq.fasta \
    > ${OUTD}/${SAMPLEID}.polished.all.hq.mapped.sam \
    2>${OUTD}/${SAMPLEID}.minimap2.log

date
