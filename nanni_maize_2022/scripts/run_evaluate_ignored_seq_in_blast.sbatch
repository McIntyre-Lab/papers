#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=eval_blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/eval_blast_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Plot distribution of cluster lengths and count how many clusters are
##     <500nt, >= 500nt and <1000nt, and >=1000nt

## Of the clusters >=1000, count how many have blast results that match
##     to a subject the “Zea mays” in the subject name

module load python/3.8

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
ID=$PROJ/compare_b73_2_mo17/blast_low_identity_PBreads
COV=$PROJ/compare_b73_2_mo17/blast_low_coverage_PBreads
OUTD=$PROJ/compare_b73_2_mo17/blast_low_ident_cov_eval
    mkdir -p ${OUTD}

#for REF in b73 mo17_yan mo17_cau; do
for REF in b73; do
    ## Remove summary counts file if it already exists
    if [[ -e ${OUTD}/summary_counts_2_${REF}.csv ]]; then
        rm ${OUTD}/summary_counts_2_${REF}.csv
    fi

    for TYPE in ident cov; do
        if [[ ${TYPE} == "ident" ]]; then
            IND=${ID}
        else
            IND=${COV}
        fi

        for FILE in $(ls ${IND}/*_2_${REF}_vs_nt.tsv); do
            SAMPLE=$(basename ${FILE} _2_${REF}_vs_nt.tsv | sed s/'blast_low_'${TYPE}'_'/''/g)
            echo "Evaluate low ${TYPE} cluster from ${SAMPLE} mapped to ${REF}"

            ## Plot distribution of sequence lengths
            python $PROJ/scripts/fasta_seq_length_distribution.py \
                -fa $PROJ/compare_b73_2_mo17/uncollapsed_${REF}_IsoSeq_reads/${SAMPLE}_unCollapsed_low_${TYPE}_${REF}_ref.fa \
                -o ${OUTD}/${SAMPLE}_unCollapsed_low_${TYPE}_${REF}_ref

            ## Count clusters and evaluate blast hits
            python $PROJ/scripts/evaluate_seq_in_blast.py \
                -l ${OUTD}/${SAMPLE}_unCollapsed_low_${TYPE}_${REF}_ref_seq_length.csv \
                -b ${FILE} \
                -f "Zea mays" \
                -f "family" \
                -f "kinase" \
                -o ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}

            ## Count the number of best hits with either "family" or "kinase" in the subject name
            COUNTFK=$(awk -F "," 'NR!=1 && ($(NF-1)==1 || $NF==1)' ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_seq_ge_1000_best_hit.csv | wc -l)

            ## Add to summary counts
            if [[ ! -e ${OUTD}/summary_counts_2_${REF}.csv ]]; then
                awk 'NR==1{print "ignored_type,sample,"$0",num_seq_ge_1000_best_hit_family_or_kinase"}' \
                    ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_counts.csv \
                    > ${OUTD}/summary_counts_2_${REF}.csv
            fi
            awk -v type=${TYPE} -v sample=${SAMPLE} -v countfk=${COUNTFK} \
                'NR!=1{print type","sample","$0","countfk}' \
                ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_counts.csv \
                >> ${OUTD}/summary_counts_2_${REF}.csv
        done
    done

done
