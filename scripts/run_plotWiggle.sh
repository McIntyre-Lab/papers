#!/bin/sh

# Set up major variables
PROJ=$MCLAB/cegs_sem_sd_paper
PYGEN=$HOME/devel/python.git/wigglePlot.py
ORIG=/mnt/storage/cegs_aln/combined
DESIGN_FILE=$PROJ/design_file/cegsV_line_list.csv

#for GENE in dsx Sxl fl\(2\)d tra tra2 Spf45 snf vir Yp1 Yp2 Yp3
for GENE in dsx inr
do
    # Fix fl2d name
    NAME=`echo $GENE | sed 's/(/_/g' | sed 's/)/_/g'`

    WORK=$HOME/tmp/cegsV_wiggle/with_fudge/$NAME
    if [[ ! -e $WORK ]]; then mkdir -p $WORK; fi;

    # One gene is too fast, so I am going to run 25 genes in each job
    for I in `seq 1 75`
    do
        # Import design information
        LINE=$(sed -n "${I}p" $DESIGN_FILE)


        python $PYGEN --gff /home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff \
                      --bam $ORIG/${LINE}*.bam \
                      -g $GENE \
                      --fudge \
                      -o $WORK/${LINE}_${NAME}.png
    done
done
echo "[`date`] Script complete"
