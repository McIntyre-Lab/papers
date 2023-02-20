#!/bin/sh


## Generate venn diagrams for counts of each feature type
## 4 groups each: Female K4, Female K27, Male K4, Male K27

module purge
module load python3/3.6
module load R/3.6

date

## Set output directory
VENN=${DABG}/venn_diagrams
    mkdir -p ${VENN}

## Generate venn diagrams
for SPECIES in ${SPECIES1} ${SPECIES2}; do

    ## Get input count file
    COUNTS=${DABG}/${SPECIES}_counts.csv

    NAMES="Female_K4,Female_K27,Male_K4,Male_K27"

    for FEATURE in 5UTR 3UTR TSS300bpWindow fragment intron intergenic; do
        VALUES=$(awk -F "," -v feature=${FEATURE} '{ \
            if($1==feature"1000"){ynnn=$2;} \
            if($1==feature"0100"){nynn=$2;} \
            if($1==feature"1100"){yynn=$2;} \
            if($1==feature"0010"){nnyn=$2;} \
            if($1==feature"1010"){ynyn=$2;} \
            if($1==feature"0110"){nyyn=$2;} \
            if($1==feature"1110"){yyyn=$2;} \
            if($1==feature"0001"){nnny=$2;} \
            if($1==feature"1001"){ynny=$2;} \
            if($1==feature"0101"){nyny=$2;} \
            if($1==feature"1101"){yyny=$2;} \
            if($1==feature"0011"){nnyy=$2;} \
            if($1==feature"1011"){ynyy=$2;} \
            if($1==feature"0001"){nnny=$2;} \
            if($1==feature"1111"){yyyy=$2;} \
            }END{print ynnn","nynn","yynn","nnyn","ynyn","nyyn","yyyn","nnny \
                ","ynny","nyny","yyny","nnyy","ynyy","nnny","yyyy}' ${COUNTS})

        python3 ${PROJ}/scripts/venn_diagram.py \
            -v ${VALUES} \
            -l ${NAMES} \
            -p ${PROJ}/scripts \
            -o ${VENN}/${SPECIES}_${FEATURE}_venn
    done

done
