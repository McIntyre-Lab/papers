#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Correct the gene_id values of the GTF output from SQANTI QC
##     to match the gene_id values associated with the transcripts
##     in the classification file

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/sqanti_post_filter_b73

python $PROJ/scripts/correct_SQANTI_QC_GTF_gene_id_02avn.py \
    -c ${IND}/SQANTI_classification.txt \
    -g ${IND}/sqanti_b73_filtered_corrected.gtf \
    -o ${IND}/sqanti_b73_filtered_corrected_associated_gene.gtf \
    > ${IND}/sqanti_b73_filtered_corrected_associated_gene_output_log.out

scp ${IND}/sqanti_b73_filtered_corrected_associated_gene.gtf \
    adalena.nanni@hpg.rc.ufl.edu:/blue/mcintyre/share/maize_ainsworth/htseq_quantify_genes_loci/.

