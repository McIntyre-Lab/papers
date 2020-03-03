#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --qos=mcintyre
#SBATCH --job-name=protein
#SBATCH --output=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB/scripts/SLURM_LOGS/protein_association_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=1

module load java/1.8.0_31 python perl

### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
PACBIO=/ufrc/mcintyre/share/maize_ainsworth/sqanti_post_filter
OUT=$PROJ/protein_db_association
    mkdir -p ${OUT}

SPECIES=maize
FULL=Zea_mays
TAXID=4577
SWISS=$PROJ/references/uniprot_Zea_mays_Swiss-Prot.fasta
TRMBL=$PROJ/references/uniprot_Zea_mays_TrEMBL.fasta
REFSEQFA=$PROJ/references/Zea_mays.refseq_proteins.faa
ENSEMBLFA=$PROJ/references/Zea_mays.B73_RefGen_v4.pep.all.fa
REFSEQGFF=$PROJ/references/ref_B73_RefGen_v4_top_level.gff3
ENSEMBLGTF=$PROJ/references/Zea_mays.B73_RefGen_v4.45.gtf
CHR=$PROJ/references/chr_accessions_B73_RefGen_v4

###############    
#### INPUT ####
###############

PACBIOGTF=${PACBIO}/sqanti_filtered_corrected.gtf
PACBIOFAA=${PACBIO}/sqanti_filtered_corrected.faa


################
#### OUTPUT ####
################

## Combined Protein Output
PROTOUT=$OUT/proteinGenomic_${SPECIES}_allRefSeqANDensembl_fast.txt

## Corrected CDS Output
CDS=$OUT/total_${SPECIES}_corrected_CDS.gtf

## Nonsense Mediated Decay Output
NMD=$OUT/total_${SPECIES}_corrected_NMD.txt

## Protein Association Output
PROTASSO=$OUT/total_${SPECIES}_corrected_proteinAssociation.txt


echo ${SWISS}
echo ${TRMBL}
echo ${ENSEMBLFA}
echo ${FULL}
echo ${TAXID}
ehco ${REFSEQGFF}
echo ${ENSEMBLGTF}
echo ${CHR}
echo ${PROTOUT}
echo ${PACBIOGTF}
echo ${PACBIOFAA}
echo ${CDS}
echo ${PROTASSO}
echo ${NMD}


## Run script to create the associations
## Script has been modified to take in the file names
SCRIPT=/ufrc/mcintyre_conesa_transvar/isoAnnot/sorghum_maize_analysis/scripts
python ${SCRIPT}/refseqEnsemblProteins2_2_mod.py \
    ${SWISS} ${TRMBL} ${ENSEMBLFA} ${FULL} ${TAXID} ${REFSEQGFF} \
    ${ENSEMBLGTF} ${CHR} ${PROTOUT} ${PACBIOGTF} ${PACBIOFAA} ${CDS} \
    ${PROTASSO} ${NMD}

