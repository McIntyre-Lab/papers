#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=comb_htseq
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/combine_htseq_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Combine individual expression files into expression matrices
##     for expected_count and TPM (isoform and gene level)

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INDIR=$PROJ/check_isoseq3_lost_genes/htseq_gene_counts
ROZ=${INDIR}/roz_combine_htseq
    mkdir -p ${ROZ}

## Merging htseq output files
##     each file has unnames columns gene_id (col 1) and count column (col 2)
##     rename count column to the sample name
##     paste all files together

for SAMP in ${INDIR}/*_ccs_gene_htseq_count.tsv; do
    NAME=$(basename ${SAMP})
    echo ${NAME}
    s=$(echo ${NAME} | sed s/_ccs_gene_htseq_count.tsv//)
    echo ${s}
    awk -v sample=${s} '{if(NR==1){print "gene_id\t"sample"\n"$1"\t"$2}else{print $1"\t"$2}}' \
        ${SAMP} > ${ROZ}/${s}.tsv
done

## paste all output together and keep gene_id and count for each sample
paste ${ROZ}/*.tsv > ${ROZ}/all.tsv

## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.csv
COLS=$(awk '{ cols = "1"; \
    for (x=1; x<=NF; x++) \
    {if (x % 2 ==0) \
        { cols = cols","x; \
        } } } END {print cols}' $ROZ/all.tsv)
echo ${COLS}

cut -f ${COLS} ${ROZ}/all.tsv > ${INDIR}/all_ccs_gene_htseq_count.tsv

rm -r ${ROZ}
