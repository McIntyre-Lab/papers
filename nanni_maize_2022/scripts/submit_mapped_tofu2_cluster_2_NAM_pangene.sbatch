#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=clust_NAM
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/cluster_NAM_%A.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1

module load python/3.6

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
REF=/blue/mcintyre/share/references/maize_b73
IND=$PROJ/compare_b73_2_mo17/sqanti_b73
OUTD=$PROJ/compare_b73_2_mo17/cluster_NAM_counts
    mkdir -p ${OUTD}

## Concatenate all classification files of mapped clusters (after tofu collapse)
## Add columns for sample and genotype
COMB=${OUTD}/all_sample_2_b73_classification.txt
head -1 ${IND}/b73/amb/*_classification.txt | \
    awk '{print "sample\tgenotype\t"$0}' \
    > ${COMB}
for FILE in $(ls ${IND}/*/*/*_classification.txt); do
    SAMPLE=$(basename ${FILE} _classification.txt)
    GENO=$(echo ${FILE} | awk -F "/" '{print $(NF-2)}')
    awk -v sample=${SAMPLE} -v geno=${GENO} \
        'NR!=1{print sample"\t"geno"\t"$0}' ${FILE} \
        >> ${COMB}
done

python $PROJ/scripts/pacbio/SQANTI_QC_assoc_gene_2_NAM_pangene_counts.py \
    -c ${COMB} \
    -v ${REF}/maizeGDB_pan_genes/B73v4_to_B73v5.tsv \
    -p ${REF}/hufford_2021_NAM/pan_gene_matrix_v3_cyverse.csv \
    -o ${OUTD}/mapped_cluster_2_NAM_pangene_counts.txt
