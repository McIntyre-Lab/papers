#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=comb_htseq
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/combine_htseq_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Combine individual expression files into expression matrices
##     for expected_count and TPM (isoform and gene level)

module purge
module load python/3.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INDIR=$PROJ/htseq_quantify_genes_loci/htseq_output_counts
ROZ=${INDIR}/roz_combine_htseq
    mkdir -p ${ROZ}

## Merging htseq output files
##     each file has unnames columns gene_id (col 1) and count column (col 2)
##     rename count column to the sample name
##     paste all files together
## run script that combines SE and PE expression counts

for i in SE PE; do
    ## expected_count
    for SAMP in ${INDIR}/*_${i}_htseq_count.tsv; do
        NAME=$(basename ${SAMP})
        echo ${NAME}
        s=$(echo ${NAME} | sed s/_${i}_htseq_count.tsv//)
        echo ${s}
        awk -v sample=${s} '{if(NR==1){print "gene_id\t"sample"\n"$1"\t"$2}else{print $1"\t"$2}}' \
            ${SAMP} > ${ROZ}/${s}_${i}.tsv
    done
    ## paste all output together and keep transcriptID and expected count for each sample
    paste ${ROZ}/*_${i}.tsv > ${ROZ}/${i}_all.tsv

    ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.csv
    COLS=$(awk '{ cols = "1"; \
        for (x=1; x<=NF; x++) \
        {if (x % 2 ==0) \
            { cols = cols","x; \
            } } } END {print cols}' $ROZ/${i}_all.tsv)
    echo ${COLS}

    cut -f ${COLS} ${ROZ}/${i}_all.tsv > ${INDIR}/${i}_htseq_count.tsv
done

## run script that combines (sums) SE and PE expected counts
python $PROJ/scripts/merge_htseq_counts.py \
    -p ${INDIR}/PE_htseq_count.tsv \
    -s ${INDIR}/SE_htseq_count.tsv \
    -o ${INDIR}/combined_htseq_count.tsv

rm ${INDIR}/PE_htseq_count.tsv
rm ${INDIR}/SE_htseq_count.tsv

rm -r ${ROZ}
