#!/bin/bash


## Generate summary file of the number of unique features
SUMMARY=${FEATuniq}/${SPECIES}_all_feature_count_summary.csv
echo "feature,unique" > ${SUMMARY}
for feature in 5UTR 3UTR TSS1kbWindow TSS300bpWindow fragments fusions introns intergenic; do
    if [[ ${feature} == "fragments" ]]; then
        BED=${FRAGMENT}
    elif [[ ${feature} == "fusions" ]]; then
        BED=${FUSION}
    elif [[ ${feature} == "introns" ]]; then
        BED=${INTRON}
    elif [[ ${feature} == "intergenic" ]]; then
        BED=${FEATuniq}/${SPECIES}_intergenic.bed
    else
        BED=${FEATuniq}/${SPECIES}_${feature}_unique.bed
    fi
    type=${feature}
    echo "${type},$(wc -l ${BED} | awk '{print $1}')" >> ${SUMMARY}
done

## Get list of unique annotations of interest for Dros PB project
echo "featureID,chrom,start,end,FBgn,FBtr,featureType,num_FBgn,flag_multigene,num_FBtr,annotation_frequency" > ${ROZ}/${SPECIES}_intergenic_unique.csv
awk '{print $4","$1","$2","$3",,,intergenic,,,,"}' ${FEATuniq}/${SPECIES}_intergenic.bed >> ${ROZ}/${SPECIES}_intergenic_unique.csv
echo "featureID,chrom,start,end,FBgn,FBtr,featureType,num_FBgn,flag_multigene,num_FBtr,annotation_frequency" > ${ROZ}/${SPECIES}_fragment_unique.csv
awk -F "," 'NR!=1{print $1","$2","$3","$4","$12","$10",fragment,"$13","$15","$11","$14}' ${FRAGMENTannot} >> ${ROZ}/${SPECIES}_fragment_unique.csv
echo "featureID,chrom,start,end,FBgn,FBtr,featureType,num_FBgn,flag_multigene,num_FBtr,annotation_frequency" > ${ROZ}/${SPECIES}_fusion_unique.csv
awk -F "," 'NR!=1{print $1","$2","$3","$4","$11","$9",fusion,"$10","$13","$8","$12}' ${FUSIONannot} >> ${ROZ}/${SPECIES}_fusion_unique.csv
echo "featureID,chrom,start,end,FBgn,FBtr,featureType,num_FBgn,flag_multigene,num_FBtr,annotation_frequency" > ${ROZ}/${SPECIES}_intron_unique.csv
awk -F "," 'NR!=1{print $1","$2","$3","$4","$6",,intron,,,,,"}' ${INTRONannot} >> ${ROZ}/${SPECIES}_intron_unique.csv

awk 'FNR==1 && NR!=1{next;}{print}' ${FEATuniq}/${SPECIES}_5UTR_unique.csv ${FEATuniq}/${SPECIES}_3UTR_unique.csv \
    ${FEATuniq}/${SPECIES}_TSS300bpWindow_unique.csv ${ROZ}/${SPECIES}_fragment_unique.csv \
    ${ROZ}/${SPECIES}_intron_unique.csv ${ROZ}/${SPECIES}_intergenic_unique.csv \
    > ${FEATuniq}/${SPECIES}_5U_3U_TSS_frag_intr_inter_uniq.csv

awk 'FNR==1 && NR!=1{next;}{print}' ${FEATuniq}/${SPECIES}_5UTR_unique.csv ${FEATuniq}/${SPECIES}_3UTR_unique.csv \
    ${FEATuniq}/${SPECIES}_TSS300bpWindow_unique.csv ${ROZ}/${SPECIES}_fusion_unique.csv \
    ${ROZ}/${SPECIES}_intron_unique.csv ${ROZ}/${SPECIES}_intergenic_unique.csv \
    > ${FEATuniq}/${SPECIES}_5U_3U_TSS_fus_intr_inter_uniq.csv
