#!/bin/bash

#### Set Directories
    PROJ=/bio/mcintyre/cegs
    WORK=$PROJ/mpileup_fb551_transcriptome

#### Because I am using an PBS Array I am pulling LINE:MV:REP from an external CSV with all possible combinations
for PBS_ARRAYID in {1..324}
do
    DESIGN_FILE=$PROJ/design_files/CEGS_57_lines_no_tech.txt
    DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)

    IFS=',' read -ra ARRAY <<< "$DESIGN"

    LINE=${ARRAY[0]}
    MV=${ARRAY[1]}
    REP=${ARRAY[2]}

    if [ ! -e $WORK/${LINE}_${MV}${REP}.mpileup ]
    then
        echo ${LINE}_${MV}${REP}
        echo $PBS_ARRAYID
    fi

done
