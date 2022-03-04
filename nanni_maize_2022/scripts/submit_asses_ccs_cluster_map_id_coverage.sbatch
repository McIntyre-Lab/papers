#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH	--job-name=count_id_cov
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/count_id_cov_%A.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=1

module load samtools/1.12

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
CLUST=$PROJ/compare_b73_2_mo17/mapping_minimap2_b73
CCS=$PROJ/check_isoseq3_lost_genes/minimap2
OUTD=$PROJ/compare_b73_2_mo17/ccs_cluster_map_id_cov
    mkdir -p ${OUTD}

CCSOUT=${OUTD}/ccs_map_id_cov_counts.csv
CLUSTOUT=${OUTD}/cluster_map_id_cov_counts.csv

## Add header to ccs and cluster count output files
header="genotype,treatment,sample,num_mapped_reads,num_map_id_ge_95,num_map_id_ge_90,num_map_id_ge_85,num_map_id_ge_80,num_map_cov_ge_95,num_map_cov_ge_90,num_map_cov_ge_85,num_map_cov_ge_80,prop_map_id_ge_95,prop_map_id_ge_90,prop_map_id_ge_85,prop_map_id_ge_80,prop_map_cov_ge_95,prop_map_cov_ge_90,prop_map_cov_ge_85,prop_map_cov_ge_80"

echo ${header} > ${CCSOUT}
echo ${header} > ${CLUSTOUT}

for FILE in $(ls ${CCS}/*_minimap2_b73.sam | grep -v unfiltered); do
    NAME=$(basename ${FILE} _minimap2_b73.sam)
    GENO=$(echo ${NAME} | awk -F "_" '{print $2}')
    TRT=$(echo ${NAME} | awk -F "_" '{print $3}')
    ## Get CCS counts
    ## NOTE: CCS sam file already filtered for supp alignments and unmapped
    ## -e expression used to filter alignments for certain values
    ##     [NM] is from the NM tag which is a count of mismatches to the reference
    ##     rlen is the length of the reference alignment
    ##     qlen is the length of the read
    CCSMAP=$(samtools view ${FILE} | wc -l)
    CCSID95=$(samtools view -e "[NM] <= (rlen*0.05)" ${FILE} | wc -l)
    CCSID90=$(samtools view -e "[NM] <= (rlen*0.10)" ${FILE} | wc -l)
    CCSID85=$(samtools view -e "[NM] <= (rlen*0.15)" ${FILE} | wc -l)
    CCSID80=$(samtools view -e "[NM] <= (rlen*0.20)" ${FILE} | wc -l)
    CCSCOV95=$(samtools view -e "(rlen / qlen) >= 0.95" ${FILE} | wc -l)
    CCSCOV90=$(samtools view -e "(rlen / qlen) >= 0.90" ${FILE} | wc -l)
    CCSCOV85=$(samtools view -e "(rlen / qlen) >= 0.85" ${FILE} | wc -l)
    CCSCOV80=$(samtools view -e "(rlen / qlen) >= 0.80" ${FILE} | wc -l)
    CCSID95PROP=$(echo ${CCSID95} ${CCSMAP} | awk '{print $1/$2}')
    CCSID90PROP=$(echo ${CCSID90} ${CCSMAP} | awk '{print $1/$2}')
    CCSID85PROP=$(echo ${CCSID85} ${CCSMAP} | awk '{print $1/$2}')
    CCSID80PROP=$(echo ${CCSID80} ${CCSMAP} | awk '{print $1/$2}')
    CCSCOV95PROP=$(echo ${CCSCOV95} ${CCSMAP} | awk '{print $1/$2}')
    CCSCOV90PROP=$(echo ${CCSCOV90} ${CCSMAP} | awk '{print $1/$2}')
    CCSCOV85PROP=$(echo ${CCSCOV85} ${CCSMAP} | awk '{print $1/$2}')
    CCSCOV80PROP=$(echo ${CCSCOV80} ${CCSMAP} | awk '{print $1/$2}')

    echo "${GENO},${TRT},${NAME},${CCSMAP},${CCSID95},${CCSID90},${CCSID85},${CCSID80},${CCSCOV95},${CCSCOV90},${CCSCOV85},${CCSCOV80},${CCSID95PROP},${CCSID90PROP},${CCSID85PROP},${CCSID80PROP},${CCSCOV95PROP},${CCSCOV90PROP},${CCSCOV85PROP},${CCSCOV80PROP}" \
        >> ${CCSOUT}

    ## Get cluster counts
    ## NOTE: cluster sam file needs filtered for supp alignments and unmapped (-F 2052)
    ## -e expression used to filter alignments for certain values
    ##     [NM] is from the NM tag which is a count of mismatches to the reference
    ##     rlen is the length of the reference alignment
    CLUSTFILE=${CLUST}/${GENO}/${TRT}/${NAME}.polished.all.hq.mapped.sam
    CLUSTMAP=$(samtools view -F 2052 ${CLUSTFILE} | wc -l)
    CLUSTID95=$(samtools view -F 2052 -e "[NM] <= (rlen*0.05)" ${CLUSTFILE} | wc -l)
    CLUSTID90=$(samtools view -F 2052 -e "[NM] <= (rlen*0.10)" ${CLUSTFILE} | wc -l)
    CLUSTID85=$(samtools view -F 2052 -e "[NM] <= (rlen*0.15)" ${CLUSTFILE} | wc -l)
    CLUSTID80=$(samtools view -F 2052 -e "[NM] <= (rlen*0.20)" ${CLUSTFILE} | wc -l)
    CLUSTCOV95=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.95" ${CLUSTFILE} | wc -l)
    CLUSTCOV90=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.90" ${CLUSTFILE} | wc -l)
    CLUSTCOV85=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.85" ${CLUSTFILE} | wc -l)
    CLUSTCOV80=$(samtools view -F 2052 -e "(rlen / qlen) >= 0.80" ${CLUSTFILE} | wc -l)
    CLUSTID95PROP=$(echo ${CLUSTID95} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTID90PROP=$(echo ${CLUSTID90} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTID85PROP=$(echo ${CLUSTID85} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTID80PROP=$(echo ${CLUSTID80} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTCOV95PROP=$(echo ${CLUSTCOV95} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTCOV90PROP=$(echo ${CLUSTCOV90} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTCOV85PROP=$(echo ${CLUSTCOV85} ${CLUSTMAP} | awk '{print $1/$2}')
    CLUSTCOV80PROP=$(echo ${CLUSTCOV80} ${CLUSTMAP} | awk '{print $1/$2}')

    echo "${GENO},${TRT},${NAME},${CLUSTMAP},${CLUSTID95},${CLUSTID90},${CLUSTID85},${CLUSTID80},${CLUSTCOV95},${CLUSTCOV90},${CLUSTCOV85},${CLUSTCOV80},${CLUSTID95PROP},${CLUSTID90PROP},${CLUSTID85PROP},${CLUSTID80PROP},${CLUSTCOV95PROP},${CLUSTCOV90PROP},${CLUSTCOV85PROP},${CLUSTCOV80PROP}" \
        >> ${CLUSTOUT}

done

