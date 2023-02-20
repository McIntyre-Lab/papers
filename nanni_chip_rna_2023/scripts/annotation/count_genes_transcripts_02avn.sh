#!/bin/bash

## Count number of genes, transcripts, fragments (unique exonic regions), fusions (blocks of exonic regions),
##     and introns in each species

SUMMgeneCount=${FEATuniq}/gene_count_summary.csv

echo "species,genes,transcripts,fragments,fusions,introns,intergenic" > ${SUMMgeneCount}

## Species1
echo "${SPECIES1},$(awk '$3=="gene"' ${GTF1} | wc -l),$(cut -f 9 ${GTF1} \
    | awk '{if($6!=""){print $6}}' | sort | uniq | wc -l),$(wc -l ${EA1}_exon_fragments_coverage.bed \
    | awk '{print $1}'),$(wc -l ${EA1}_fusions_coverage.bed | awk '{print $1}'),$(wc -l ${EA1}_introns_from_fusions.bed \
    | awk '{print $1}'),$(wc -l ${FEATuniq}/${SPECIES1}_intergenic.bed | awk '{print $1}')" \
    >> ${SUMMgeneCount}

## Species2
echo "${SPECIES2},$(awk '$3=="gene"' ${GTF2} | wc -l),$(cut -f 9 ${GTF2} \
    | awk '{if($6!=""){print $6}}' | sort | uniq | wc -l),$(wc -l ${EA2}_exon_fragments_coverage.bed \
    | awk '{print $1}'),$(wc -l ${EA2}_fusions_coverage.bed | awk '{print $1}'),$(wc -l ${EA2}_introns_from_fusions.bed \
    | awk '{print $1}'),$(wc -l ${FEATuniq}/${SPECIES2}_intergenic.bed | awk '{print $1}')" \
    >> ${SUMMgeneCount}
