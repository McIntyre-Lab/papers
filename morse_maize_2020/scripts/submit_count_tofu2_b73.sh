#!/bin/bash

### counts reads from tofu2/cupcake2 output

## Check argument number
if [ "$#" -ne 1 ]; then
    echo "Must give one argument of tofu2 file path including prefix.
	For example : /ufrc/mcintyre_conesa_transvar/dros_PacBio/tofu2/mel/mel
		for mel.cluster_report.csv and mel.collapsed.group.txt in that path"
    exit 1
fi

## Check argument is a file
if [[ -f "${1}.cluster_report.csv" && -f "${1}.collapsed.group.txt" ]]; then
    PREFIX=${1}
    echo " prefix is " $PREFIX
    NAME=`basename "$PREFIX"`
    echo " basename is " $NAME

### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    OUT=$PROJ/pacbio_analysis/counting_recheck
        if [ ! -e $OUT ]; then mkdir -p $OUT; fi

    ## Count reads contained in clusters via output from isoseq3_make_cluster_report.py
    wc -l ${PREFIX}.cluster_report.csv | awk -v name=$NAME '{ print name  "\t" $0-1}'  > $OUT/${NAME}_tofu2_b73_counts_cluster_report.txt

    ## Count unique isoforms from collapse_isoforms_by_sam.py
    wc -l ${PREFIX}.collapsed.group.txt | awk -v name=$NAME '{print name "\t" $1}' > $OUT/${NAME}_tofu2_b73_counts_collaped_group.txt

else
    echo "Must give one argument of tofu2 file path including prefix.
        For example : /ufrc/mcintyre_conesa_transvar/dros_PacBio/tofu2/mel/mel
                for mel.cluster_report.csv and mel.collapsed.group.txt in that path"
    exit 1
fi

