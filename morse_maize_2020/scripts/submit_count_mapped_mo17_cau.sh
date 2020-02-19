#!/bin/bash

### counts reads from alignment output (regardless of program used as long as it is a bam or sam file)

module load samtools/1.3.1

## Check argument number
if [ "$#" -ne 1 ]; then
    echo "Must give one argument of pacbio alignment file path to bam or sam output.
	For example : Proj/isoseq3_analysis/mapping/mel/mel_combined.polished.all.hq.mapped.sam"
    exit 1
fi

## Check argument is a file
if [[ -f $1  ]]; then
    SAM=${1}
else
    echo "Must give one argument of pacbio alignment file path to bam or sam output.
        For example : Proj/isoseq3_analysis/mapping/mel/mel_combined.polished.all.hq.mapped.sam"
    exit 1
fi

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUT=$PROJ/pacbio_analysis/counting_recheck
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi

samtools view -F 2052 $SAM | wc -l | awk -OFS="\t" -v sam=$SAM '{print sam "\t" $0}' > $OUT/$(basename ${SAM})_minimap_mo17_cau_qc_cnt.csv

