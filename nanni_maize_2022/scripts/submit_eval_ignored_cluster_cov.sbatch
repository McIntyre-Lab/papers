#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH	--job-name=eval_clust_cov
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/eval_clust_cov_%A.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=4

module load samtools/1.12 python/3.6

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/compare_b73_2_mo17/ccs_cluster_map_id_cov
    mkdir -p ${OUTD}

for REF in b73 mo17_yan mo17_cau; do
    ## Get input directories
    IND=$PROJ/compare_b73_2_mo17/mapping_minimap2_${REF}
    IGNORE=$PROJ/compare_b73_2_mo17/uncollapsed_${REF}_IsoSeq_reads
    IGNORETYPE="low_cov"
    ROZ=$OUTD/roz_ignored_2_${REF}
        mkdir -p ${ROZ}

    CLUSTOUT=${OUTD}/ignored_cluster_map_cov_counts_2_${REF}.csv

    ## Add header to ccs and cluster count output files
    header="genotype,treatment,sample,ignore_type,num_reads,num_map_cov_lt_10,num_map_cov_lt_15,num_map_cov_lt_20,num_map_cov_lt_25,num_map_cov_lt_30,num_map_cov_lt_35,num_map_cov_lt_40,num_map_cov_lt_45,num_map_cov_lt_50,num_map_cov_lt_55,num_map_cov_lt_60,num_map_cov_lt_65,num_map_cov_lt_70,num_map_cov_lt_75,prop_map_cov_lt_10,prop_map_cov_lt_15,prop_map_cov_lt_20,prop_map_cov_lt_25,prop_map_cov_lt_30,prop_map_cov_lt_35,prop_map_cov_lt_40,prop_map_cov_lt_45,prop_map_cov_lt_50,prop_map_cov_lt_55,prop_map_cov_lt_60,prop_map_cov_lt_65,prop_map_cov_lt_70,prop_map_cov_lt_75"

    echo ${header} > ${CLUSTOUT}

    for FULLFILE in $(ls ${IND}/*/*/*.polished.all.hq.mapped.sam); do

        NAME=$(basename ${FULLFILE} .polished.all.hq.mapped.sam)
        GENO=$(echo ${NAME} | awk -F "_" '{print $2}')
        TRT=$(echo ${NAME} | awk -F "_" '{print $3}')

        ## Subset SAM file for ignored cluster IDs
        python $PROJ/scripts/Common_extracFeatureFromList.py \
            -l ${IGNORE}/list_${NAME}.ignored_low_cov_ids.${REF}_ref.txt \
            -i ${FULLFILE} \
            -t sam \
            -o ${ROZ}/${NAME}.ignored_low_cov.sam
        FILE=${ROZ}/${NAME}.ignored_low_cov.sam

        ## Get cluster counts
        ## NOTE: cluster sam file needs filtered for supp alignments and unmapped (-F 2052)
        ## -e expression used to filter alignments for certain values
        ##     [NM] is from the NM tag which is a count of mismatches to the reference
        ##     rlen is the length of the reference alignment
        TOTAL=$(samtools view -F 2052 ${FILE} | wc -l)
        CLUSTCOV10=$(samtools view -F 2052 -e "(rlen / qlen) < 0.10" ${FILE} | wc -l)
        CLUSTCOV15=$(samtools view -F 2052 -e "(rlen / qlen) < 0.15" ${FILE} | wc -l)
        CLUSTCOV20=$(samtools view -F 2052 -e "(rlen / qlen) < 0.20" ${FILE} | wc -l)
        CLUSTCOV25=$(samtools view -F 2052 -e "(rlen / qlen) < 0.25" ${FILE} | wc -l)
        CLUSTCOV30=$(samtools view -F 2052 -e "(rlen / qlen) < 0.30" ${FILE} | wc -l)
        CLUSTCOV35=$(samtools view -F 2052 -e "(rlen / qlen) < 0.35" ${FILE} | wc -l)
        CLUSTCOV40=$(samtools view -F 2052 -e "(rlen / qlen) < 0.40" ${FILE} | wc -l)
        CLUSTCOV45=$(samtools view -F 2052 -e "(rlen / qlen) < 0.45" ${FILE} | wc -l)
        CLUSTCOV50=$(samtools view -F 2052 -e "(rlen / qlen) < 0.50" ${FILE} | wc -l)
        CLUSTCOV55=$(samtools view -F 2052 -e "(rlen / qlen) < 0.55" ${FILE} | wc -l)
        CLUSTCOV60=$(samtools view -F 2052 -e "(rlen / qlen) < 0.60" ${FILE} | wc -l)
        CLUSTCOV65=$(samtools view -F 2052 -e "(rlen / qlen) < 0.65" ${FILE} | wc -l)
        CLUSTCOV70=$(samtools view -F 2052 -e "(rlen / qlen) < 0.70" ${FILE} | wc -l)
        CLUSTCOV75=$(samtools view -F 2052 -e "(rlen / qlen) < 0.75" ${FILE} | wc -l)
        CLUSTCOV10PROP=$(echo ${CLUSTCOV10} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV15PROP=$(echo ${CLUSTCOV15} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV20PROP=$(echo ${CLUSTCOV20} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV25PROP=$(echo ${CLUSTCOV25} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV30PROP=$(echo ${CLUSTCOV30} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV35PROP=$(echo ${CLUSTCOV35} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV40PROP=$(echo ${CLUSTCOV40} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV45PROP=$(echo ${CLUSTCOV45} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV50PROP=$(echo ${CLUSTCOV50} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV55PROP=$(echo ${CLUSTCOV55} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV60PROP=$(echo ${CLUSTCOV60} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV65PROP=$(echo ${CLUSTCOV65} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV70PROP=$(echo ${CLUSTCOV70} ${TOTAL} | awk '{print $1/$2}')
        CLUSTCOV75PROP=$(echo ${CLUSTCOV75} ${TOTAL} | awk '{print $1/$2}')

        echo "${GENO},${TRT},${NAME},${IGNORETYPE},${TOTAL},${CLUSTCOV10},${CLUSTCOV15},${CLUSTCOV20},${CLUSTCOV25},${CLUSTCOV30},${CLUSTCOV35},${CLUSTCOV40},${CLUSTCOV45},${CLUSTCOV50},${CLUSTCOV55},${CLUSTCOV60},${CLUSTCOV65},${CLUSTCOV70},${CLUSTCOV75},${CLUSTCOV10PROP},${CLUSTCOV15PROP},${CLUSTCOV20PROP},${CLUSTCOV25PROP},${CLUSTCOV30PROP},${CLUSTCOV35PROP},${CLUSTCOV40PROP},${CLUSTCOV45PROP},${CLUSTCOV50PROP},${CLUSTCOV55PROP},${CLUSTCOV60PROP},${CLUSTCOV65PROP},${CLUSTCOV70PROP},${CLUSTCOV75PROP}" \
            >> ${CLUSTOUT}
    done
    rm -r ${ROZ}
done

