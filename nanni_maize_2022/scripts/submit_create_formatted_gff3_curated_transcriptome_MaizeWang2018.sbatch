#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=format_gff
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/create_formatted_gff3_MaizeWang2018_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb               # per processor memory
#SBATCH --array=1-5

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### Shell script to convert GTF/GFF to Flybase-style GFF3

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

## Set input and output directories
    INDIR=$PROJ/MaizeWang2018_transcriptome_eval/subset_FSM_ISM_NIC_NNC
    OUTDIR=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome
    mkdir -p ${OUTDIR}
    echo -e "GTF/GFF formatting for ${SAMPLEID}...\nOutput in ${OUTDIR}"

### Path to GTF/GFF file to convert
    GFF=${INDIR}/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.gtf

### Output GFF path and name. Make sure this ends with ".gff"
    GFFOUT=${OUTDIR}/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.converted.gff

### Path to Event Analysis install. This is the folder containing "docs", "src", run_buildAnnotations.sh, etc.
    EVENTDIR=$PROJ/scripts/event_analysis
    SCRIPTS=${EVENTDIR}/src

### Set temp directory
    TMPDIR=${OUTDIR}/roz_${SAMPLEID}_gtf
    if [[ ! -e ${TMPDIR} ]]; then mkdir -p ${TMPDIR}; fi

#### Checking if a gff.db file exists, and if not then create one
    echo "Checking if user-supplied GTF/GFF file has a pre-generated database file"
        if [ ! -e ${GFF}.db ]
        then
        echo "No database file found! Generating..."
        python $SCRIPTS/make_gff_db.py --gff $GFF
        echo "Database generated!"
        else
        echo "Database file found! Skipping generation."
        fi

### Convert user-supplied GTF/GFF file into GFF3 format
    echo "Converting user-supplied GTF/GFF to GFF3 format"

    if [ ! -e ${GFFOUT} ]
    then
    python $SCRIPTS/convertGTF2GFF3.py --input $GFF --output ${GFFOUT}
    sort -k1,1 -k4n -k5n ${GFFOUT} > ${TMPDIR}/temp.gff
    mv ${TMPDIR}/temp.gff ${GFFOUT}
    else
    echo "Converted file found! Skipping conversion."
    fi

    if [ ! -e ${GFFOUT}.db ]
    then
    python $SCRIPTS/make_gff_db.py --gff ${GFFOUT}
    else
    echo "Converted database file found! Skipping generation."
    fi

    rm -r ${TMPDIR}

