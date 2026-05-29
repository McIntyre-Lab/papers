#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=trand2bed.%j
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/trand2bed.%j.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/trand2bed.%j.err
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b


### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=${PROJ}/trand_gene_isoseq_FSM_ISM_NIC_02ndk
OUTD=${IND}
ROZ=${OUTD}/roz
    mkdir -p ${ROZ}

##Cat all fragment and region EA files together
##Frags
(head -n 1 ${IND}/1_isoseq_FSM_ISM_NIC_event_analysis_ef.csv && tail -n +2 -q ${IND}/*_isoseq_FSM_ISM_NIC_event_analysis_ef.csv) > ${ROZ}/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv


#Regions
(head -n 1 ${IND}/1_isoseq_FSM_ISM_NIC_event_analysis_er.csv && tail -n +2 -q ${IND}/*_isoseq_FSM_ISM_NIC_event_analysis_er.csv) > ${ROZ}/allChr_isoseq_FSM_ISM_NIC_event_analysis_er.csv


##Convert TranD EA fragments output to bed file for quantificatio

#FRAGMENTS

##Trand ouput columns are:
##gene_id,er_id,ef_id,ef_chr,ef_start,ef_end,ef_strand,ef_exon_ids,ef_transcript_ids,exons_per_ef,transcripts_per_ef,ef_ir_flag,ea_annotation_frequency

##For bed file
##ef_chr   ef_start   ef_end   ef_id
awk -F "," 'NR>1 {print $4"\t"$5"\t"$6"\t"$3}'  ${ROZ}/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv > ${OUTD}/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.bed

##FUSIONS

#TranD output columsn
##gene_id,er_id,er_chr,er_start,er_end,er_strand,er_exon_ids,er_transcript_ids,gene_transcript_ids,exons_per_er,transcripts_per_er,transcripts_per_gene,er_ir_flag,er_annotation_frequency

##For bed file
##er_chr   er_start   er_end   er_id

awk -F "," 'NR>1 {print $3"\t"$4"\t"$5"\t"$2}' ${ROZ}/allChr_isoseq_FSM_ISM_NIC_event_analysis_er.csv > ${OUTD}/allChr_isoseq_FSM_ISM_NIC_event_analysis_er.bed
