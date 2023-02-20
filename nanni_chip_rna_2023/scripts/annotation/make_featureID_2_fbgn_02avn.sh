#!/bin/bash

## Generate featureID to FBgn file

COMBINED=${FEATuniq}/${SPECIES}_featureID_2_fbgn.csv
echo "featureID,FBgn" > ${COMBINED}
for FEAT in 5UTR 3UTR TSS1kbWindow fragment fusion intron; do

    if [[ ${FEAT} == "fragment" ]]; then
        ## Print fragment IDs and associated FBgn
        cut -d "," -f 1,12 ${FRAGMENTannot} \
            | awk -F "," '{if($1!="fragment_id"){split($2,a,"|"); \
            for(i=1;i<=length(a);i++){print $1","a[i]}}}' | sort | uniq \
            >> ${COMBINED}
    elif [[ ${FEAT} == "fusion" ]]; then
        ## Print fusion IDs and associated FBgn
        cut -d "," -f 1,11 ${FUSIONannot} \
            | awk -F "," '{if($1!="fragment_id"){split($2,a,"|"); \
            for(i=1;i<=length(a);i++){print $1","a[i]}}}' | sort | uniq \
            >> ${COMBINED}
    elif [[ ${FEAT} == "intron" ]]; then
        ## Print intron IDs and associated FBgn
        cut -d "," -f 1,6 ${INTRONannot} \
       	    | awk -F "," '{if($1!="intron_id"){split($2,a,"|"); \
       	    for(i=1;i<=length(a);i++){print $1","a[i]}}}' | sort | uniq \
       	    >> ${COMBINED}
    else
        INPUT=${FEATuniq}/${SPECIES}_${FEAT}_xcrpt_annotation.csv
        ## Print feature IDs and associated FBgn
        awk -F "," 'NR!=1{print $7","$4}' ${INPUT} \
            | sort | uniq >> ${COMBINED}
    fi

done
