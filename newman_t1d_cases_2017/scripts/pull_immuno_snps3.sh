#!/bin/bash

# Script to pull Immunochip SNPs for each Immunobase gene

## Set files

#PROJ=/home/jrbnewman/concannon/eqtl_analysis
PROJ=/mnt/data/eqtls


VCF=$PROJ/original_data/eurfam_vcf_pruned.vcf
BED=$PROJ/immunogenes.bed
LIST=$PROJ/immunogenes.txt
OUT=$PROJ/immunobase_snps


## Iterate through gene list

for GENE in $(cat $LIST);
    do
    echo "#necessary_header" > temp.bed
    awk -v gene="$GENE" '{if ($4 == gene) print $0 }' $BED >> temp.bed
    vcftools --vcf $VCF --out $OUT/${GENE} --bed temp.bed --recode
    tail -n+30 $OUT/${GENE}.recode.vcf | awk -v gene="$GENE" '{print gene "\t" $0}' > $OUT/$GENE.vcf
    rm $OUT/${GENE}.recode.vcf
done

rm temp.bed
rm $OUT/*.log



