#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH	--job-name=count_id_cov
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/count_id_cov_%A.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=4

module load samtools/1.12

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
CLUST=$PROJ/compare_b73_2_mo17/mapping_minimap2_b73
CCS=$PROJ/check_isoseq3_lost_genes/minimap2
OUTD=$PROJ/compare_b73_2_mo17/ccs_cluster_map_id_cov
    mkdir -p ${OUTD}

for REF in b73 mo17_yan mo17_cau; do
    ## Get input directories
#    CLUST=$PROJ/compare_b73_2_mo17/mapping_minimap2_${REF}
    CCS=$PROJ/check_isoseq3_lost_genes/minimap2

    CCSOUT=${OUTD}/ccs_map_id_cov_counts_2_${REF}.csv
#    CLUSTOUT=${OUTD}/cluster_map_id_cov_counts_2_${REF}.csv

    ## Add header to ccs and cluster count output files
    header="genotype,treatment,sample,num_mapped_reads,num_map_id_ge_95,num_map_id_ge_90,num_map_id_ge_85,num_map_id_ge_80,num_map_cov_ge_95,num_map_cov_ge_90,num_map_cov_ge_85,num_map_cov_ge_80,num_map_id_cov_95,prop_map_id_ge_95,prop_map_id_ge_90,prop_map_id_ge_85,prop_map_id_ge_80,prop_map_cov_ge_95,prop_map_cov_ge_90,prop_map_cov_ge_85,prop_map_cov_ge_80,prop_map_id_cov_95,num_unmapped_reads"

    echo ${header} > ${CCSOUT}
#    echo ${header} > ${CLUSTOUT}

    for FILE in $(ls ${CCS}/*_minimap2_${REF}.sam | grep -v unfiltered); do
        NAME=$(basename ${FILE} _minimap2_${REF}.sam)
        GENO=$(echo ${NAME} | awk -F "_" '{print $2}')
        TRT=$(echo ${NAME} | awk -F "_" '{print $3}')
        ## Get CCS counts
        ## NOTE: CCS sam file already filtered for supp alignments and unmapped
        ## -e expression used to filter alignments for certain values
        ##     [NM] is from the NM tag which is a count of mismatches to the reference
        ##     rlen is the length of the reference alignment
        ##     qlen is the length of the read
        CCSUNMAP=$(samtools view -f 4 ${CCS}/${NAME}_unfiltered_minimap2_${REF}.sam | wc -l)
        CCSMAP=$(samtools view ${FILE} | wc -l)
        CCSID95=$(samtools view -e "[NM] <= (rlen*0.05)" ${FILE} | wc -l)
        CCSID90=$(samtools view -e "[NM] <= (rlen*0.10)" ${FILE} | wc -l)
        CCSID85=$(samtools view -e "[NM] <= (rlen*0.15)" ${FILE} | wc -l)
        CCSID80=$(samtools view -e "[NM] <= (rlen*0.20)" ${FILE} | wc -l)
        CCSCOV95=$(samtools view -e "(rlen / qlen) >= 0.95" ${FILE} | wc -l)
        CCSCOV90=$(samtools view -e "(rlen / qlen) >= 0.90" ${FILE} | wc -l)
        CCSCOV85=$(samtools view -e "(rlen / qlen) >= 0.85" ${FILE} | wc -l)
        CCSCOV80=$(samtools view -e "(rlen / qlen) >= 0.80" ${FILE} | wc -l)
        CCSIDCOV95=$(samtools view -e "[NM] <= (rlen*0.05) && (rlen / qlen) >= 0.95" ${FILE} | wc -l)
        CCSID95PROP=$(echo ${CCSID95} ${CCSMAP} | awk '{print $1/$2}')
        CCSID90PROP=$(echo ${CCSID90} ${CCSMAP} | awk '{print $1/$2}')
        CCSID85PROP=$(echo ${CCSID85} ${CCSMAP} | awk '{print $1/$2}')
        CCSID80PROP=$(echo ${CCSID80} ${CCSMAP} | awk '{print $1/$2}')
        CCSCOV95PROP=$(echo ${CCSCOV95} ${CCSMAP} | awk '{print $1/$2}')
        CCSCOV90PROP=$(echo ${CCSCOV90} ${CCSMAP} | awk '{print $1/$2}')
        CCSCOV85PROP=$(echo ${CCSCOV85} ${CCSMAP} | awk '{print $1/$2}')
        CCSCOV80PROP=$(echo ${CCSCOV80} ${CCSMAP} | awk '{print $1/$2}')
        CCSIDCOV95PROP=$(echo ${CCSIDCOV95} ${CCSMAP} | awk '{print $1/$2}')

        echo "${GENO},${TRT},${NAME},${CCSMAP},${CCSID95},${CCSID90},${CCSID85},${CCSID80},${CCSCOV95},${CCSCOV90},${CCSCOV85},${CCSCOV80},${CCSIDCOV95},${CCSID95PROP},${CCSID90PROP},${CCSID85PROP},${CCSID80PROP},${CCSCOV95PROP},${CCSCOV90PROP},${CCSCOV85PROP},${CCSCOV80PROP},${CCSIDCOV95PROP},${CCSUNMAP}" \
            >> ${CCSOUT}

        ## Get cluster counts
        ## NOTE: cluster sam file needs filtered for supp alignments and unmapped (-F 2052)
        ## -e expression used to filter alignments for certain values
        ##     [NM] is from the NM tag which is a count of mismatches to the reference
        ##     rlen is the length of the reference alignment
#        CLUSTFILE=${CLUST}/${GENO}/${TRT}/${NAME}.polished.all.hq.mapped.sam
#        CLUSTUNMAP=$(awk '$2==4' ${CLUST}/${GENO}/${TRT}/${NAME}.polished.all.hq.mapped.sorted.sam | wc -l)
#        CLUSTMAP=$(samtools view -F 2052 ${CLUSTFILE} | wc -l)
#        CLUSTID95=$(samtools view -F 2052 -e "[NM] <= (rlen*0.05)" ${CLUSTFILE} | wc -l)
#        CLUSTID90=$(samtools view -F 2052 -e "[NM] <= (rlen*0.10)" ${CLUSTFILE} | wc -l)
#        CLUSTID85=$(samtools view -F 2052 -e "[NM] <= (rlen*0.15)" ${CLUSTFILE} | wc -l)
#        CLUSTID80=$(samtools view -F 2052 -e "[NM] <= (rlen*0.20)" ${CLUSTFILE} | wc -l)
#        CLUSTCOV95=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.95" ${CLUSTFILE} | wc -l)
#        CLUSTCOV90=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.90" ${CLUSTFILE} | wc -l)
#        CLUSTCOV85=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.85" ${CLUSTFILE} | wc -l)
#        CLUSTCOV80=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.80" ${CLUSTFILE} | wc -l)
#        CLUSTIDCOV95=$(samtools view -F 2052 -e "[NM] <= (rlen*0.05) && (rlen / qlen) >= 0.95" ${CLUSTFILE} | wc -l)
#        CLUSTID95PROP=$(echo ${CLUSTID95} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTID90PROP=$(echo ${CLUSTID90} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTID85PROP=$(echo ${CLUSTID85} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTID80PROP=$(echo ${CLUSTID80} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTCOV95PROP=$(echo ${CLUSTCOV95} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTCOV90PROP=$(echo ${CLUSTCOV90} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTCOV85PROP=$(echo ${CLUSTCOV85} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTCOV80PROP=$(echo ${CLUSTCOV80} ${CLUSTMAP} | awk '{print $1/$2}')
#        CLUSTIDCOV95PROP=$(echo ${CLUSTIDCOV95} ${CLUSTMAP} | awk '{print $1/$2}')

#        echo "${GENO},${TRT},${NAME},${CLUSTMAP},${CLUSTID95},${CLUSTID90},${CLUSTID85},${CLUSTID80},${CLUSTCOV95},${CLUSTCOV90},${CLUSTCOV85},${CLUSTCOV80},${CLUSTIDCOV95},${CLUSTID95PROP},${CLUSTID90PROP},${CLUSTID85PROP},${CLUSTID80PROP},${CLUSTCOV95PROP},${CLUSTCOV90PROP},${CLUSTCOV85PROP},${CLUSTCOV80PROP},${CLUSTIDCOV95PROP},${CLUSTUNMAP}" \
#            >> ${CLUSTOUT}
    done
done

## To output in the format of STable 1 do the following
## CCS
#echo "

#STable 1 % of mapped CCS with >=95% id (top) and >=95% coverage (bottom) to B73 / Mo17 YAN / Mo17 CAU reference

#"
#for SAMPLE in $(awk -F "," 'NR!=1{print $3}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv); do
#    echo -e "\n\n${SAMPLE}"
#    ## column 3 = sample
#    ## column 14 = prop_map_id_ge_95
#    ## column 18 = prop_map_cov_ge_95
#    ## column 23 = num_unmapped_reads
#    ## column 4 = num_mapped_reads
#    B73MAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
#    YANMAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_yan.csv)
#    CAUMAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_cau.csv)
#    B73ID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
#    CAUID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_cau.csv)
#    YANID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_yan.csv)
#    B73COV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
#    CAUCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_cau.csv)
#    YANCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/ccs_map_id_cov_counts_2_mo17_yan.csv)
#    echo "${B73MAP} ${YANMAP} ${CAUMAP}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#    echo "${B73ID} ${YANID} ${CAUID}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#    echo "${B73COV} ${YANCOV} ${CAUCOV}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#done


#echo "

#STable 1 % of mapped clusters with >=95% id OR >=95% coverage to B73 / Mo17 YAN / Mo17 CAU reference

#"
#for SAMPLE in $(awk -F "," 'NR!=1{print $3}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv); do
#    echo -e "\n\n${SAMPLE}"
#    ## column 3 = sample
#    ## column 14 = prop_map_id_ge_95
#    ## column 18 = prop_map_cov_ge_95
#    ## column 23 = num_unmapped_reads
#    ## column 4 = num_mapped_reads
#    B73MAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv)
#    YANMAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_yan.csv)
#    CAUMAP=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $4/($23 + $4)}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_cau.csv)
#    B73ID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv)
#    CAUID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_cau.csv)
#    YANID=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $14}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_yan.csv)
#    B73COV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv)
#    CAUCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_cau.csv)
#    YANCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $18}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_yan.csv)
#    echo "${B73MAP} ${YANMAP} ${CAUMAP}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#    echo "${B73ID} ${YANID} ${CAUID}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#    echo "${B73COV} ${YANCOV} ${CAUCOV}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
#done


echo "

STable 1 % of mapped CCS with >=95% id AND >=95% coverage to B73 / Mo17 YAN / Mo17 CAU reference

"
for SAMPLE in $(awk -F "," 'NR!=1{print $3}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv); do
    echo -e "\n\n${SAMPLE}"
    ## column 3 = sample
    ## column 22 = prop_map_id_cov_95
    B73IDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
    CAUIDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
    YANIDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/ccs_map_id_cov_counts_2_b73.csv)
    echo "${B73IDCOV} ${YANIDCOV} ${CAUIDCOV}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
done


echo "

STable 1 % of mapped clusters with >=95% id AND >=95% coverage to B73 / Mo17 YAN / Mo17 CAU reference

"
for SAMPLE in $(awk -F "," 'NR!=1{print $3}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv); do
    echo -e "\n\n${SAMPLE}"
    ## column 3 = sample
    ## column 22 = prop_map_id_cov_95
    B73IDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/cluster_map_id_cov_counts_2_b73.csv)
    CAUIDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_cau.csv)
    YANIDCOV=$(awk -F "," -v sample=${SAMPLE} '$3==sample{print $22}' ${OUTD}/cluster_map_id_cov_counts_2_mo17_yan.csv)
    echo "${B73IDCOV} ${YANIDCOV} ${CAUIDCOV}" | awk '{printf("%.2f% \/ %.2f% \/ %.2f%\n", $1 * 100, $2 * 100, $3 * 100)}' 2>/dev/null
done

