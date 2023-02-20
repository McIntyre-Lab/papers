#!/bin/bash


## Identify regions outside of the genic features
## including 5'UTR, 3'UTR, TSS (1kb region cntered on 5'UTR start),
##     fragments (unique exonic regions), and introns

## Get bed file of full genome (use chrom.sizes)
awk '{print $1"\t"0"\t"$2}' ${CHROM} > ${ROZ}/${SPECIES}.chrom.sizes.bed

## Combine coordinates of unique features (5'UTR,3'UTR,TSS1kbWindow,TSS300bpWindow,fragments,fusions,intronic)
cat ${FEATuniq}/*_unique.bed ${FRAGMENT} ${FUSION} ${INTRON} | \
    sort -k1,1 -k2,2n > ${FEATuniq}/${SPECIES}_all_features.bed

## Get intergenic regions using bedtools
bedtools subtract \
    -a ${ROZ}/${SPECIES}.chrom.sizes.bed \
    -b ${FEATuniq}/${SPECIES}_all_features.bed \
    > ${ROZ}/temp_${SPECIES}_intergenic.bed

## Add intergenic featureIDs
## Require a region length >50bp - <50bp are places in another file in roz
awk '$3-$2>50{print $1"\t"$2"\t"$3"\tintergenic_"$1"_"$2"_"$3}' \
    ${ROZ}/temp_${SPECIES}_intergenic.bed \
    > ${FEATuniq}/${SPECIES}_intergenic.bed
awk '$3-$2<=50{print $1"\t"$2"\t"$3"\tintergenic_"$1"_"$2"_"$3}' \
    ${ROZ}/temp_${SPECIES}_intergenic.bed \
    > ${ROZ}/${SPECIES}_intergenic_LE_50.bed

## Remove temprary files
rm ${ROZ}/temp_${SPECIES}_intergenic.bed
rm ${ROZ}/${SPECIES}.chrom.sizes.bed
