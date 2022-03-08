#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=eval_blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/eval_blast_2_clusters_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Of the clusters >=1000, count how many have blast result to clusters

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
    if [[ -e ${OUTD}/summary_counts_2_${REF}_vs_b73_mo17_clusters.csv ]]; then
        rm ${OUTD}/summary_counts_2_${REF}_vs_b73_mo17_clusters.csv
    fi

    for TYPE in ident cov; do
        if [[ ${TYPE} == "ident" ]]; then
            IND=${ID}
        else
            IND=${COV}
        fi

        for COMPARE in b73 mo17; do
            for FILE in $(ls ${IND}/*_2_${REF}_vs_${COMPARE}_collapsed_clusters.tsv); do
                SAMPLE=$(basename ${FILE} _2_${REF}_vs_${COMPARE}_collapsed_clusters.tsv | sed s/'blast_low_'${TYPE}'_'/''/g)
                echo "Evaluate low ${TYPE} cluster from ${SAMPLE} mapped to ${REF} to ${COMPARE} collapsed clusters"

                ## Count clusters and evaluate blast hits
                python $PROJ/scripts/evaluate_seq_in_blast.py \
                    -l ${OUTD}/${SAMPLE}_unCollapsed_low_${TYPE}_${REF}_ref_seq_length.csv \
                    -b ${FILE} \
                    -o ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_2_${COMPARE}

                ## Add to summary counts
                if [[ ! -e ${OUTD}/summary_counts_2_${REF}_vs_b73_mo17_clusters.csv ]]; then
                    awk 'NR==1{print "ignored_type,compared_clusters,sample,"$0}' \
                        ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_2_${COMPARE}_counts.csv \
                        > ${OUTD}/summary_counts_2_${REF}_vs_b73_mo17_clusters.csv
                fi
                awk -v type=${TYPE} -v compare=${COMPARE} -v sample=${SAMPLE} \
                    'NR!=1{print type","compare","sample","$0}' \
                    ${OUTD}/${SAMPLE}_2_${REF}_low_${TYPE}_2_${COMPARE}_counts.csv \
                    >> ${OUTD}/summary_counts_2_${REF}_vs_b73_mo17_clusters.csv
            done
        done
    done
done
