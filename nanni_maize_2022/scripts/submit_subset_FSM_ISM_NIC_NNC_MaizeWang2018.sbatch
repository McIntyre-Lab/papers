#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=subset
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/subset_FSM_ISM_NIC_NNC_MaizeWang2018_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

## Identify and remove unspliced fragments and
##     select representative reference transcripts for FSM/ISM
## 1) Flag monoexon reads in classification file
## 2) Select NIC/NNC that are not unspliced fragments
## 3) Subset SQANTI GTF (with ref corrected gene_id) for NIC/NNC that pass selection
## 4) Get unique FSM/ISM associated reference transcripts
## 5) Subset Reference GTF for FSM/ISM associated transcripts
## 6) Concatenate GTF files of (3) and (5)


module load python/3.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
SCRIPTS=$PROJ/scripts/pacbio_transcriptome_eval
OUTD=$PROJ/MaizeWang2018_transcriptome_eval
SQANTIOUT=${OUTD}/sqanti_qc
SUBSETOUT=${OUTD}/subset_FSM_ISM_NIC_NNC
    mkdir -p ${SUBSETOUT}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Set log file
LOG=${SUBSETOUT}/${SAMPLEID}_subset_counts.log
ROZ=${SUBSETOUT}/roz_${SAMPLEID}_subset
    mkdir -p ${ROZ}

## Get reference genome FASTA and GTF
if [[ ${SPECIES} == "maize" ]]; then
    REF_GTF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
else
    REF_GTF=/blue/mcintyre/share/references/sorghum_bicolor/ensembl_NCBIv3.49/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.49.gtf
fi


date

echo "
*** Remove Unspliced NIC/NNC Fragments ***

Inputs: ${SQANTIOUT}/${SAMPLEID}_classification.txt
        ${SQANTIOUT}/${SAMPLEID}_corrected_associated_gene.gtf
"
## Flag monoexon reads in classification file
python ${SCRIPTS}/flag_monoexon_02avn.py \
    -c ${SQANTIOUT}/${SAMPLEID}_classification.txt \
    -o ${SUBSETOUT}/${SAMPLEID}_corrected_associated_gene_monoexon_freq.txt \
    > ${LOG}
echo -e "\t...${SAMPLEID} monoexon reads flagged"

## Select NIC/NNC that are not unspliced fragments
##     Monoexon PB to multiexon ref (will be excluded)
##     Multiexon PB to multiexon ref (keep)
##     Monoexon PB to monoexon ref (keep)
## All monoexon to monoexon are FSM or ISM, there are no multiexon to monoexon
## $NF = flag_ref_monoexon, $(NF-1) = flag_PB_monoexon
awk '($6=="novel_in_catalog" || $6=="novel_not_in_catalog") && \
    $NF==0 && $(NF-1)==0{print $1}' \
    ${SQANTIOUT}/${SAMPLEID}_classification_flag_mono.txt \
    > ${ROZ}/temp_pbID.txt

## Subset SQANTI GTF (with ref corrected gene_id) for NIC/NNC that pass selection
python ${SCRIPTS}/subset_gtf_simple_02avn.py \
    -g ${SQANTIOUT}/${SAMPLEID}_corrected_associated_gene.gtf \
    -t transcript_id \
    -i ${ROZ}/temp_pbID.txt \
    -o ${ROZ}/sqanti_nic_nnc.gtf

## Get unique FSM/ISM associated reference transcripts
awk '$6=="full-splice_match" || $6=="incomplete-splice_match"{print $8}' \
    ${SQANTIOUT}/${SAMPLEID}_classification_flag_mono.txt \
    | sort | uniq > ${ROZ}/temp_ref.txt

## Subset Reference GTF for FSM/ISM associated transcripts
python ${SCRIPTS}/subset_gtf_simple_02avn.py \
    -g ${REF_GTF} \
    -t transcript_id \
    -i ${ROZ}/temp_ref.txt \
    -o ${ROZ}/ref_fsm_ism.gtf

## Concatenate GTF files of FSM/ISM associated reference transcripts
##     and NIC/NNC selected transcripts
cat ${ROZ}/ref_fsm_ism.gtf ${ROZ}/sqanti_nic_nnc.gtf \
    > ${SUBSETOUT}/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.gtf

## Count
echo "
$(awk 'NR!=1' ${SQANTIOUT}/${SAMPLEID}_classification_flag_mono.txt | wc -l) total reads
$(cat ${ROZ}/temp_ref.txt ${ROZ}/temp_pbID.txt | wc -l) after subsetting
     $(cat ${ROZ}/temp_ref.txt | wc -l) FSM/ISM associated reference transcripts
     $(cat ${ROZ}/temp_pbID.txt | wc -l) NIC/NNC that pass removal of unspliced fragments
">${LOG}

rm -r ${ROZ}

date
