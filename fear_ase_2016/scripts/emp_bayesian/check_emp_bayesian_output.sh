#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: Compare the number of lines that were analyzable with the
#   number of lines that were analyzed and make sure these numbers matched.
# 
#===============================================================================

PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
INPUT=$PROJ/emp_bayesian/input/ase_dataset_for_bayesian.csv
OUTPUT=$PROJ/emp_bayesian/output/PG_emp_bayesian_results.csv

# Logs to output counts
MYLOG=$PROJ/emp_bayesian/output/PG_emp_bayesian_results_size_check.log

# Count how many fusions were going to be analyzed
INCOUNT=`awk 'BEGIN{FS=","; count=0}
             {
                 if( $16 == 1 )
                 { count+=1 }
             } 
             END{print count}' $INPUT`

    if [ -e $OUTPUT ]
    then
        OUTCOUNT=`wc -l $OUTPUT | cut -f1 -d' '`
        NOHEAD=$(expr $OUTCOUNT - 1)

        if [[ $INCOUNT == $NOHEAD ]]
        then
            printf "OK:\t$INCOUNT\t$NOHEAD\n" >> $MYLOG
        else
            printf "BAD COUNT:\t$INCOUNT\t$NOHEAD\n" >> $MYLOG
        fi
    else
        printf "NO OUTPUT: $OUTPUT\n" >> $MYLOG
    fi
