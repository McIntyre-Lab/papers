#!/bin/bash

## Generate featureID to featureType file

COMBINED=${FEATuniq}/${SPECIES}_all_featureID_2_featureType.csv
echo "featureID,featureType" > ${COMBINED}
for FEAT in 5UTR 3UTR TSS1kbWindow TSS300bpWindow fragment fusion intron intergenic; do
    if [[ ${FEAT} == "fragment" ]]; then
        if [[ ${SPECIES} == ${SPECIES1} ]]; then
            BED=${EA1}_exon_fragments_coverage.bed
        else
            BED=${EA2}_exon_fragments_coverage.bed
        fi
    elif [[ ${FEAT} == "fusion" ]]; then
        if [[ ${SPECIES} == ${SPECIES1} ]]; then
            BED=${EA1}_fusions_coverage.bed
        else
            BED=${EA2}_fusions_coverage.bed
        fi
    elif [[ ${FEAT} == "intron" ]]; then
        if [[ ${SPECIES} ==	${SPECIES1} ]]; then
            BED=${EA1}_introns_from_fusions.bed
        else
            BED=${EA2}_introns_from_fusions.bed
       	fi
    elif [[ ${FEAT} == "intergenic" ]]; then
        BED=${FEATuniq}/${SPECIES}_${FEAT}.bed
    else
        BED=${FEATuniq}/${SPECIES}_${FEAT}_unique.bed
    fi
    
    awk -v feature=${FEAT} '{print $4","feature}' ${BED} >> ${COMBINED}
done
