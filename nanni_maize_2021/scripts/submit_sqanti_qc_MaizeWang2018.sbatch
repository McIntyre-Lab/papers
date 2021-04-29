#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/sqanti_qc_MaizeWang2018_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

### SQANTI isoseq QC of maize PacBio data -
### 9 samples total

module purge
module load pacbio/5.1.0 isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval
TOFUOUT=${OUTD}/tofu2
SQANTIOUT=${OUTD}/sqanti_qc
    mkdir -p ${SQANTIOUT}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Get reference genome FASTA and GTF
if [[ ${SPECIES} == "maize" ]]; then
    REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
    REF_GTF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
else
    REF=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
    REF_GTF=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.49.gtf
fi

date

## GTF (is GTF but has .gff ending) from output of collapse_isoforms_by_sam.py
PACBIO_GTF=${TOFUOUT}/${SAMPLEID}.collapsed.gff
ABUNDANCE=${TOFUOUT}/${SAMPLEID}.collapsed.abundance.txt

date
echo "***SQANTI***
"
python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py \
    -g ${PACBIO_GTF} ${REF_GTF} \
    ${REF} -o ${SAMPLEID} \
    -d ${SQANTIOUT} -fl ${ABUNDANCE} 2>${SQANTIOUT}/${SAMPLEID}.sqanti.log

date
