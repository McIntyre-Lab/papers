#!/bin/bash

## Count SQANTI QC classifications in CCS reads from each sample
##     mapped to the 3 references

## Set directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/compare_b73_2_mo17

## Add header to counts file
COUNTS=${IND}/sqanti_ccs_read_counts_all_refs.csv
echo "ref,sample,num_FSM,num_ISM,num_NIC,num_NNC,total_read" \
    > ${COUNTS}

## Loop over references to count SQANTI QC outputs
for REF in b73 mo17_yan mo17_cau; do
#for REF in b73 mo17_yan; do
    SQANTI=${IND}/sqanti_${REF}_ccs_reads
    ## Loop over samples
    for SAMPLE in $(ls ${SQANTI}); do
        ## Count number of reads in each category
        ##     structural_category in column 6
        CLASS=${SQANTI}/${SAMPLE}/${SAMPLE}_ccs_classification.txt
        FSM=$(awk '$6=="full-splice_match"' ${CLASS} | wc -l)
        ISM=$(awk '$6=="incomplete-splice_match"' ${CLASS} | wc -l)
        NIC=$(awk '$6=="novel_in_catalog"' ${CLASS} | wc -l)
        NNC=$(awk '$6=="novel_not_in_catalog"' ${CLASS} | wc -l)
        TOTAL=$(awk 'NR!=1' ${CLASS} | wc -l)

        ## Output to counts file
        echo "${REF},${SAMPLE},${FSM},${ISM},${NIC},${NNC},${TOTAL}" \
            >> ${COUNTS}
    done
done


## To output in the format of STable 1 do the following
#for SAMPLE in $(ls ${SQANTI}); do
#    echo -e "\n\n${SAMPLE}"
#    awk -F "," -v sample=${SAMPLE} \
#        '$2==sample{if($1=="b73"){FSM1=$3; ISM1=$4; NIC1=$5; NNC1=$6}\
#        if($1=="mo17_yan"){FSM2=$3; ISM2=$4; NIC2=$5; NNC2=$6}\
#        if($1=="mo17_cau"){FSM3=$3; ISM3=$4; NIC3=$5; NNC3=$6}}\
#        END{print FSM1" / "FSM2" / "FSM3"\n"ISM1" / "ISM2" / "ISM3"\n"NNC1" / "NNC2" / "NNC3"\n"NIC1" / "NIC2" / "NIC3"\n"}' \
#        ${COUNTS}
#done
