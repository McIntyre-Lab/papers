#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=build_annot
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/buildAnnotations_MaizeWang2018_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb               # per processor memory
#SBATCH --array=1-5

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### Shell script to build Event Analysis annotations

## Set project directory
    PROJ=/blue/mcintyre/share/maize_ainsworth

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

### Path to Event Analysis install
    EVENTDIR=$PROJ/scripts/event_analysis

### Label to prefix to annotation files (e.g. mm10_refseq, dmel617, hg38_ens)
    PREFIX=${SAMPLEID}

### Output directory. If it does not exist, the annotation build script with create it
    OUTDIR=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/${SAMPLEID}_150bp_annotations

    echo -e "Event Analysis annotations for ${SAMPLEID}...\nOutput in ${OUTDIR}"

### Path to formatted GFF3 file (should in a FlyBase GFF3 format)
    GFF=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.converted.gff

### Path to genome FASTA
if [[ ${SPECIES} == "maize" ]]; then
    FASTA=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
else
    FASTA=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
fi
### Size (in nt) of reads or each read pair (set to maximum if you have reads of
### various lengths, i.e. from adapter trimming)
    READSIZE=150

### Build Event Analysis annotations
    cd ${EVENTDIR}
    sh ./run_buildAnnotations.sh ${PREFIX} ${GFF} ${FASTA} ${READSIZE} ${OUTDIR}

### Remove the fasta files of EA output to reduce memory
rm ${OUTDIR}/*.fa
