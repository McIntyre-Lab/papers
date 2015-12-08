#!/bin/bash

# Script to pull Immunochip SNPs for each Immunobase gene

## Set files

PROJ=/home/jrbnewman/concannon/eqtl_analysis

VCF=$PROJ/original_data/eurfam_vcf2.vcf
BED=$PROJ/scripts/immunobase_genes.bed
LIST=$PROJ/immunogenes.txt
OUT=$PROJ/immunobase3


## Iterate through gene list

for GENE in $(cat $LIST);
    do
    echo "#necessary_header" > temp.bed
    grep $GENE $BED >> temp.bed
    vcftools --vcf $VCF --out $OUT/${GENE} --bed temp.bed --recode
    tail -n+30 $OUT/${GENE}.recode.vcf | awk -v gene="$GENE" '{print gene "\t" $0}' > $OUT/$GENE.vcf
    rm $OUT/${GENE}.recode.vcf
done

rm *.bed
rm $OUT/*.log



