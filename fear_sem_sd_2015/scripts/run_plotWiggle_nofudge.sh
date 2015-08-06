#!/bin/sh

# Set up major variables
PROJ=$MCLAB/cegs_sem_sd_paper
PYGEN=$HOME/devel/python.git/wigglePlot.py
ORIG=/mnt/storage/cegs_aln/combined
DESIGN_FILE=$PROJ/design_file/cegsV_line_list.csv

#for GENE in InR dsx Sxl fl\(2\)d tra tra2 Spf45 snf vir Yp1 Yp2 Yp3 Psa
for GENE in Psa
do
    # Fix fl2d name
    NAME=`echo $GENE | sed 's/(/_/g' | sed 's/)/_/g'`

    WORK=$HOME/tmp/cegsV_wiggle/no_fudge/$NAME
    if [[ ! -e $WORK ]]; then mkdir -p $WORK; fi;

    # Iterate over all lines
    for I in `seq 1 75`
    do
        # Import design information
        LINE=$(sed -n "${I}p" $DESIGN_FILE)

        # Create a smple id that matches the VCF sample names
        SAMPLE=`echo $LINE | sed 's/r/Raleigh_/g' | sed 's/w/Winters_/g'`

        python $PYGEN --bam $ORIG/${LINE}*.bam \
                      --gff $HOME/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff \
                      --gene $GENE \
                      --vcf $HOME/storage/useful_dmel_data/CEGS.68.lines.raw.SNPs.filtered.set.1.recode.vcf.gz \
                      --sample $SAMPLE \
                      -o $WORK/${LINE}_${NAME}.png
    done
done
echo "[`date`] Script complete"
