#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=length_FSM
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/length_FSM_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Get lengths of all ccs reads (after trimming)

## Plot distribution of ccs read and reference transcript lengths

module load gffread/0.12.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/compare_b73_2_mo17/sqanti_b73_ccs_reads
OUTD=$PROJ/compare_b73_2_mo17/ccs_read_length
    mkdir -p ${OUTD}
ROZ=${OUTD}/roz_FSM_length
    mkdir -p ${ROZ}

## Get reference transcript fasta
REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
REF_GTF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf

## Get reference transcript FA from GTF and genome sequence
gffread -w ${ROZ}/ref_b73.fa \
    -g ${REF} ${REF_GTF}

## Set up env to have seaborn v0.9.0 for python code to run properly
export PATH=/blue/mcintyre/adalena.nanni/conda_envs/plot_seaborn_0.9.0/bin:$PATH

## Set up counts file
COUNTS=${OUTD}/all_vs_ref_FSM_ISM_length_counts.csv
echo "genotype,num_ccs,num_FSM_ccs,num_FSM_assoc_ref,num_FSM_ccs_5perc_of_ref,num_FSM_ccs_10perc_of_ref,num_FSM_ccs_15perc_of_ref,num_FSM_ccs_20perc_of_ref,prop_total_ccs_FSM_5perc_of_ref,prop_total_ccs_FSM_10perc_of_ref,prop_total_ccs_FSM_15perc_of_ref,prop_total_ccs_FSM_20perc_of_ref,num_ISM_ccs,num_ISM_assoc_ref,num_ISM_ccs_5perc_of_ref,num_ISM_ccs_10perc_of_ref,num_ISM_ccs_15perc_of_ref,num_ISM_ccs_20perc_of_ref,prop_total_ccs_ISM_5perc_of_ref,prop_total_ccs_ISM_10perc_of_ref,prop_total_ccs_ISM_15perc_of_ref,prop_total_ccs_FSM_20perc_of_ref" \
    > ${COUNTS}

## Concatenate samples within each genotype
for GENO in b73 c123 hp301 mo17 nc338; do
    echo "Evaluate ${GENO} FSM ccs reads"
    ## Get ccs read length from SQANTI QC classification file
    ## col 1 = isoform
    ## col 4 = length
    ## col 6 = structural_category
    echo "name,length" > ${ROZ}/${GENO}_FSM_ccs.csv
    for FILE in $(ls ${IND}/*${GENO}*/*_classification.txt); do
        awk 'NR!=1 && $6 == "full-splice_match"{print $1","$4}' ${FILE} \
            >> ${ROZ}/${GENO}_FSM_ccs.csv
    done

    ## Get associated B73 reference transcript lengths from SQANTI QC classification file
    ## col 8 = associated_transcript
    ## col 9 = ref_length
    echo "name,length" > ${ROZ}/${GENO}_FSM_assoc_b73_ref.csv
    for	FILE in	$(ls ${IND}/*${GENO}*/*_classification.txt); do
        awk 'NR!=1 && $6 == "full-splice_match"{print $8","$9}' ${FILE}
    done | sort | uniq >> ${ROZ}/${GENO}_FSM_assoc_b73_ref.csv

    ## Plot distribution of FSM sequence lengths vs. reference lengths
    echo "...plot"
    if [[ ! -e ${OUTD}/ref_b73_seq_length.csv ]]; then
        python $PROJ/scripts/plot_seq_length_distribution.py \
            -l ${ROZ}/${GENO}_FSM_assoc_b73_ref.csv \
            -fa ${ROZ}/ref_b73.fa \
            -pdf \
            -o ${OUTD} \
            -p ${GENO}_FSM_ref_vs_b73_ref
    else
        python $PROJ/scripts/plot_seq_length_distribution.py \
            -l ${ROZ}/${GENO}_FSM_assoc_b73_ref.csv \
            -l ${OUTD}/ref_b73_seq_length.csv \
            -pdf \
            -o ${OUTD} \
            -p ${GENO}_FSM_ref_vs_b73_ref
    fi

    echo "...plot"
    python $PROJ/scripts/plot_seq_length_distribution.py \
        -l ${ROZ}/${GENO}_FSM_ccs.csv \
        -o ${OUTD}

    ## Plot scatter plot of all ccs lengths vs. associated reference transcripts
    ## Get only ccs reads that have an associated reference transcript (FSM/ISM)
    echo "isoform,length,structural_category,associated_transcript,ref_length" \
        > ${ROZ}/${GENO}_vs_ref_classification.csv
    for FILE in $(ls ${IND}/*${GENO}*/*_classification.txt); do
        awk -F "\t" 'NR!=1 && $8!="novel"{print $1","$4","$6","$8","$9}' \
            ${FILE} \
            >> ${ROZ}/${GENO}_vs_ref_classification.csv
    done
    python $PROJ/scripts/plot_scatter_from_CSV.py \
        -i ${ROZ}/${GENO}_vs_ref_classification.csv \
        -x ref_length \
        -xlab "B73 Reference Transcript Length" \
        -y length \
        -ylab "${GENO} CCS Read Length" \
        -g structural_category \
        -o ${OUTD}/${GENO}_FSM_ISM_vs_ref_length.png

    ## Plot scatter plot of FSM only
    echo "isoform,length,structural_category,associated_transcript,ref_length" \
        > ${ROZ}/${GENO}_FSM_vs_ref_classification.csv
    for FILE in $(ls ${IND}/*${GENO}*/*_classification.txt); do
        awk -F "\t" 'NR!=1 && $6 == "full-splice_match"{print $1","$4","$6","$8","$9}' \
            ${FILE} \
            >> ${ROZ}/${GENO}_FSM_vs_ref_classification.csv
    done
    python $PROJ/scripts/plot_scatter_from_CSV.py \
        -i ${ROZ}/${GENO}_FSM_vs_ref_classification.csv \
        -x ref_length \
        -xlab "B73 Reference Transcript Length" \
        -y length \
        -ylab "${GENO} FSM CCS Read Length" \
        -g structural_category \
        -o ${OUTD}/${GENO}_FSM_vs_ref_length.png

    ## Get list of FSM associated reference transcripts and the number of ccs with each
    ##     sorted by length of the reference transcript (longest first)
    awk '$6=="full-splice_match"{print $8" "$9}' ${IND}/*${GENO}*/*_classification.txt | \
        sort | uniq -c | sort -k3,3nr | \
        awk 'BEGIN{print "num_ccs,assoc_transcript,len_assoc_transcript"}{print $1","$2","$3}' \
        > ${OUTD}/${GENO}_FSM_assoc_ref_length.csv

    ## Plot scatter plot of ISM only
    echo "isoform,length,structural_category,associated_transcript,ref_length" \
        > ${ROZ}/${GENO}_ISM_vs_ref_classification.csv
    for FILE in $(ls ${IND}/*${GENO}*/*_classification.txt); do
        awk -F "\t" 'NR!=1 && $6 == "incomplete-splice_match"{print $1","$4","$6","$8","$9}' \
            ${FILE} \
            >> ${ROZ}/${GENO}_ISM_vs_ref_classification.csv
    done
    python $PROJ/scripts/plot_scatter_from_CSV.py \
        -i ${ROZ}/${GENO}_ISM_vs_ref_classification.csv \
        -x ref_length \
        -xlab "B73 Reference Transcript Length" \
        -y length \
        -ylab "${GENO} ISM CCS Read Length" \
        -g structural_category \
        -o ${OUTD}/${GENO}_ISM_vs_ref_length.png

    ## Get list of ISM associated reference transcripts	and the	number of ccs with each
    ## 	   sorted by length of the reference transcript	(longest first)
    awk '$6=="incomplete-splice_match"{print $8" "$9}' ${IND}/*${GENO}*/*_classification.txt | \
        sort | uniq -c | sort -k3,3nr | \
        awk 'BEGIN{print "num_ccs,assoc_transcript,len_assoc_transcript"}{print $1","$2","$3}' \
	> ${OUTD}/${GENO}_ISM_assoc_ref_length.csv

    ## Get counts and proportion of ccs reads that are FSM/ISM and within X% of the reference length
    TOTAL=$(awk 'FNR==2 && NR!=1{next;}{print}' ${IND}/*${GENO}*/*_classification.txt | wc -l)
    FSM=$(awk 'NR!=1' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | wc -l)
    FSMref=$(awk -F "," 'NR!=1{print $4}' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | sort | uniq | wc -l)
    FSM5=$(awk -F "," '$2>=($5*0.95) && $2<=($5*1.05)' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | wc -l)
    FSM10=$(awk -F "," '$2>=($5*0.90) && $2<=($5*1.10)' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | wc -l)
    FSM15=$(awk -F "," '$2>=($5*0.85) && $2<=($5*1.15)' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | wc -l)
    FSM20=$(awk -F "," '$2>=($5*0.80) && $2<=($5*1.20)' ${ROZ}/${GENO}_FSM_vs_ref_classification.csv | wc -l)
    FSM5prop=$(echo "${FSM5} ${TOTAL}" | awk '{print $1/$2}')
    FSM10prop=$(echo "${FSM10} ${TOTAL}" | awk '{print $1/$2}')
    FSM15prop=$(echo "${FSM15} ${TOTAL}" | awk '{print $1/$2}')
    FSM20prop=$(echo "${FSM20} ${TOTAL}" | awk '{print $1/$2}')
    ISM=$(awk 'NR!=1' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | wc -l)
    ISMref=$(awk -F "," 'NR!=1{print $4}' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | sort | uniq | wc -l)
    ISM5=$(awk -F "," '$2>=($5*0.95) && $2<=($5*1.05)' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | wc -l)
    ISM10=$(awk -F "," '$2>=($5*0.90) && $2<=($5*1.10)' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | wc -l)
    ISM15=$(awk -F "," '$2>=($5*0.85) && $2<=($5*1.15)' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | wc -l)
    ISM20=$(awk -F "," '$2>=($5*0.80) && $2<=($5*1.20)' ${ROZ}/${GENO}_ISM_vs_ref_classification.csv | wc -l)
    ISM5prop=$(echo "${ISM5} ${TOTAL}" | awk '{print $1/$2}')
    ISM10prop=$(echo "${ISM10} ${TOTAL}" | awk '{print $1/$2}')
    ISM15prop=$(echo "${ISM15} ${TOTAL}" | awk '{print $1/$2}')
    ISM20prop=$(echo "${ISM20} ${TOTAL}" | awk '{print $1/$2}')

    echo "${GENO},${TOTAL},${FSM},${FSMref},${FSM5},${FSM10},${FSM15},${FSM20},${FSM5prop},${FSM10prop},${FSM15prop},${FSM20prop},${ISM},${ISMref},${ISM5},${ISM10},${ISM15},${ISM20},${ISM5prop},${ISM10prop},${ISM15prop},${ISM20prop}" \
        >> ${COUNTS}

done

rm -r ${ROZ}
