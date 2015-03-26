#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./create_wiggle_plots.sh 
# 
#   DESCRIPTION: R cannot handle the size of these data sets, so I need to
#   parse the data to the gene level and then feed that to R.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 10/16/2012 12:23:00 PM EDT
#      REVISION:  ---
#===============================================================================

# Directories
PROJ=$MCLAB/arbeitman_fru_network
INDIR=$PROJ/data/for_wiggles
OUTDIR=$PROJ/reports_internal/wiggles
TMPDIR=$HOME/tmp/wiggles

if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi
if [ ! -d $TMPDIR ]; then mkdir $TMPDIR; fi

# Files
RWIG=$PROJ/r_programs/wiggleplots_v3.R
COORD=/home/jfear/mclab/useful_dmel_data/flybase530/symbol2coord.csv

# Gene of Interest
for GENE in Sxl fru dsx tra orb2 Fbp1 Fbp2
#for GENE in Sxl
do
    # Get Gene Coordinates

    LINE=`grep "$GENE" $COORD`

    IFS=',' read -ra ARRAY <<< "$LINE"

    CHROM=${ARRAY[2]}
    START=`echo "${ARRAY[3]} - 1000" | bc`
    END=`echo "${ARRAY[4]} + 1000" | bc`

    cd $INDIR

    for SEX in MALE FEMALE
    #for SEX in MALE 
    do
        BGA='';
        BGB='';
        BGC='';

        for FRU in A B C
        do
            FILE=AH_${SEX}_FRUM_${FRU}.csv

            if [ ! -e $TMPDIR/${GENE}_AH_${SEX}_FRUM_${FRU}.csv ]
            then
                awk '/'$CHROM'/ {if('$START' <= $2 && '$END' >= $2) {print $0}}' FS=, $FILE  > $TMPDIR/${GENE}_AH_${SEX}_FRUM_${FRU}.csv
            fi

            CURRBG=BG$FRU
            
            eval ${CURRBG}=`grep $FILE $INDIR/background.csv | cut -d',' -f3`

        done

        #Rscript $RWIG $TMPDIR $GENE $SEX $BGA $BGB $BGC

        #rm $TMPDIR/*.csv
    done
done
