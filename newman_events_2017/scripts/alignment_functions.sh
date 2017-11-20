#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: A series of bowtie and last functions for use with the
#   alignment pipeline.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 07/30/2013 11:53:38 AM EDT
#      REVISION:  ---
#===============================================================================

function check_bowtie {

    echo "START Bowtie Size Checker [`date`]" >>$MYLOG
    # Do a size check and store the result in VAR
    VAR=`/scratch/lfs/mcintyre/python.git/SizeChecker.py -i $READS -aln $ALN -unaln $UNALN`

    # Check VAR; if it is -1, then bowtie produced a different output size than input size

    if [[ $VAR -eq -1 ]]
    then
        echo $NAME "- WARNING: Bowtie input/output size inconsistent \n" >>$ALN_ERROR_LOG
    fi

    echo "FINISH Bowtie Size Checker [`date`]" >>$MYLOG
}


function bowtie_se_uniq {
    cd $TMPDIR

    echo "START Single End bowtie [`date`]" >>$MYLOG
    bowtie -S \
           $btqual \
           --best \
           --tryhard \
           --strata \
           --chunkmbs 1024 \
           -p $NUMPROCS \
           -m 1 \
           -v 3 \
           --un ${NAME}_unaln_bt.fq \
           --max ${NAME}_ambig_bt.fq \
           $REF \
           $READS \
           2>$ALNLOGS/${NAME}_bowtie.log \
           | perl -ne 'if(/^@/){next} @s = split; if($s[1] != 4) {print}' \
           > $TMPDIR/${NAME}_aln_bt.sam \
           2>>$MYLOG
        echo "FINISH Single End bowtie [`date`]" >>$MYLOG

        if [ ! -e ${NAME}_unaln_bt.fq ]
        then
            echo "WARNING: There were no unaligned reads" >> $MYLOG
        fi

        ALN=${NAME}_aln_bt.sam

        UNALN=UNALN.fq
        cat ${NAME}_unaln_bt.fq ${NAME}_ambig_bt.fq > $UNALN
        
        check_bowtie

        READS=${NAME}_unaln_bt.fq
        rm $UNALN

}

function bowtie_se_all {
    cd $TMPDIR

    echo "START Single End bowtie [`date`]" >>$MYLOG
    bowtie -S \
           $btqual \
           --best \
           --tryhard \
           --strata \
           --chunkmbs 1024 \
           -p $NUMPROCS \
           -a \
           -v 3 \
           --un ${NAME}_unaln_bt.fq \
           $REF \
           $READS \
           2>$ALNLOGS/${NAME}_bowtie.log \
           | perl -ne 'if(/^@/){next} @s = split; if($s[1] != 4) {print}' \
           > ${NAME}_aln_bt.sam \
           2>>$MYLOG
        echo "FINISH Single End bowtie [`date`]" >>$MYLOG

        if [ ! -e ${NAME}_unaln_bt.fq ]
        then
            echo "WARNING: There were no unaligned reads" >> $MYLOG
        fi

        ALN=${NAME}_aln_bt.sam

        UNALN=UNALN.fq
        cat ${NAME}_unaln_bt.fq > $UNALN
        
        check_bowtie

        READS=${NAME}_unaln_bt.fq
        rm $UNALN

}

function bowtie_pe_uniq {
    cd $TMPDIR

    echo "START Paired End bowtie [`date`]" >>$MYLOG
    bowtie -S \
           $btqual \
           --best \
           --tryhard \
           --strata \
           --chunkmbs 1024 \
           -p $NUMPROCS \
           -m 1 \
           -v 3 \
           --un ${NAME}_unaln_bt.fq \
           --max ${NAME}_ambig_bt.fq \
           $REF \
           -1 $READ1 \
           -2 $READ2 \
           2>$ALNLOGS/${NAME}_bowtie.log \
           | perl -ne 'if(/^@/){next} @s = split; if($s[1] != 77 && $s[1] != 141) {print}' \
           >${NAME}_aln_bt.sam \
           2>>$MYLOG
        echo "FINISH Paired End bowtie [`date`]" >>$MYLOG

        if [ ! -e ${NAME}_unaln_bt_1.fq ] || [ ! -e ${NAME}_unaln_bt_2.fq ]
        then
            echo "WARNING: There were no unaligned reads" >> $MYLOG
        fi

        READS=READS.fq
        cat $READ1 $READ2 > $READS

        ALN=${NAME}_aln_bt.sam $ALN

        UNALN=UNALN.fq
        cat ${NAME}_unaln_bt_*.fq ${NAME}_ambig_bt_*.fq > $UNALN

        du -sh *.fq >>$MYLOG
        
        check_bowtie
        rm $UNALN $READS

        READ1=${NAME}_unaln_bt_1.fq
        READ2=${NAME}_unaln_bt_2.fq
}


function bowtie_pe_all {
    cd $TMPDIR

    echo "START Paired End bowtie [`date`]" >>$MYLOG
    bowtie -S \
           $btqual \
           --best \
           --tryhard \
           --strata \
           --chunkmbs 1024 \
           -p $NUMPROCS \
           -a  \
           -v 3 \
           --un ${NAME}_unaln_bt.fq \
           $REF \
           -1 $READ1 \
           -2 $READ2 \
           2>$ALNLOGS/${NAME}_bowtie.log \
           | perl -ne 'if(/^@/){next} @s = split; if($s[1] != 77 && $s[1] != 141) {print}' \
           >${NAME}_aln_bt.sam \
           2>>$MYLOG
        echo "FINISH Paired End bowtie [`date`]" >>$MYLOG

        if [ ! -e ${NAME}_unaln_bt_1.fq ] || [ ! -e ${NAME}_unaln_bt_2.fq ]
        then
            echo "WARNING: There were no unaligned reads" >> $MYLOG
        fi

        READS=READS.fq
        cat $READ1 $READ2 > $READS

        ALN=${NAME}_aln_bt.sam

        UNALN=UNALN.fq
        cat ${NAME}_unaln_bt_*.fq > $UNALN

        du -sh *.fq >>$MYLOG
        
        check_bowtie
        rm $UNALN $READS

        READ1=${NAME}_unaln_bt_1.fq
        READ2=${NAME}_unaln_bt_2.fq
}


function last_se_uniq {
    cd $TMPDIR

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READS $NUMPROCS
        rm $READS

    # START LAST for each of the NUMPROCS new -split files
        echo "START Single End LAST [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Single End LAST [`date`]" >>$MYLOG

        rm *-split.fq

    # Convert MAF to SAM
        echo "START Converting from MAF to SAM [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
            maf-convert.py sam ${I}_last.maf > ${I}_last.sam 2>>$MYLOG &
        done
        wait

        rm *_last.maf

        echo "FINISH Converting to SAM [`date`]" >>$MYLOG

        for ((I=0; I<$NUMPROCS; I++))
        do
            cat ${I}_last.sam >> ${NAME}_last.sam
            rm ${I}_last.sam
        done

    # PARSE LAST FILES
        echo "START Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG
        perl $PROJ/scripts/parse_last_sam_v2.pl \
            ${NAME}_last.sam \
            ${NAME}_ambig_last.sam \
            ${NAME}_uniq_last.sam \
            > $ALNLOGS/${NAME}_LAST.log

        rm ${NAME}_last.sam ${NAME}_ambig_last.sam 

        echo "FINISH Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG

}


function last_se_all {
    cd $TMPDIR

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READS $NUMPROCS
        rm $READS

    # START LAST for each of the NUMPROCS new -split files
        echo "START Single End LAST [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Single End LAST [`date`]" >>$MYLOG

        rm *-split.fq

    # Convert MAF to SAM
        echo "START Converting from MAF to SAM [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
            maf-convert.py sam ${I}_last.maf > ${I}_last.sam 2>>$MYLOG &
        done
        wait

        rm *_last.maf

        echo "FINISH Converting to SAM [`date`]" >>$MYLOG

        for ((I=0; I<$NUMPROCS; I++))
        do
            cat ${I}_last.sam >> ${NAME}_last.sam
            rm ${I}_last.sam
        done

    # PARSE LAST FILES
        echo "START Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG
        perl $PROJ/scripts/parse_last_sam_v2.pl \
            ${NAME}_last.sam \
            ${NAME}_ambig_last.sam \
            ${NAME}_uniq_last.sam \
            > $ALNLOGS/${NAME}_LAST.log

        rm ${NAME}_last.sam 

        echo "FINISH Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG

}

function last_pe_uniq {
    cd $TMPDIR

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READ1 $NUMPROCS

    # START LAST for each of the NUMPROCS new -split files
        echo "START Paired End LAST READ 1 [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Paired End LAST READ 1 [`date`]" >>$MYLOG

        cat *_last.maf > read1.maf
        rm *_last.maf *-split.fq

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READ2 $NUMPROCS

    # START LAST for each of the NUMPROCS new -split files

        echo "START Paired End LAST READ 2 [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 20 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Paired End LAST READ 2 [`date`]" >>$MYLOG

        cat *_last.maf > read2.maf
        rm *_last.maf *-split.fq

    # Combine Reads
        echo "START Combining Paired End  MAFs [`date`]" >>$MYLOG
            last-pair-probs.py -r read1.maf read2.maf > ${NAME}_last.maf 2>>$MYLOG
        echo "FINISH Combining Paired End  MAFs [`date`]" >>$MYLOG

    # Convert MAF to SAM
        echo "START Converting from MAF to SAM [`date`]" >>$MYLOG
            maf-convert.py sam ${NAME}_last.maf > ${NAME}_last.sam 2>>$MYLOG 
        echo "FINISH Converting to SAM [`date`]" >>$MYLOG

    # PARSE LAST FILES
        echo "START Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG
        perl $PROJ/scripts/parse_last_sam_v2.pl \
            ${NAME}_last.sam \
            ${NAME}_ambig_last.sam \
            ${NAME}_uniq_last.sam \
            > $ALNLOGS/${NAME}_LAST.log

        rm ${NAME}_last.sam ${NAME}_ambig_last.sam 

        echo "FINISH Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG

}


function last_pe_all {
    cd $TMPDIR

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READ1 $NUMPROCS
        rm $READ1

    # START LAST for each of the NUMPROCS new -split files
        echo "START Paired End LAST READ 1 [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Paired End LAST READ 1 [`date`]" >>$MYLOG

        cat *_last.maf > read1.maf
        rm *_last.maf *-split.fq

    # The following Python script turns READS into NUMPROCS separate files, going from 0-split.fq to (NUMPROCS-1)-split.fq. 
    # READS is the location where the output files are placed

        python /scratch/lfs/mcintyre/python.git/SeparateLast.py $READ2 $NUMPROCS
        rm $READ2

    # START LAST for each of the NUMPROCS new -split files
        echo "START Paired End LAST READ 2 [`date`]" >>$MYLOG
        for ((I=0; I<$NUMPROCS; I++))
        do
                READS=${I}-split.fq

                lastal \
                    -l 25 \
                    -Q $lastqual \
                    $LASTREF \
                    $READS \
                    2>>$MYLOG \
                    >${I}_last.maf &
        done
        wait

        echo "FINISH Paired End LAST READ 2 [`date`]" >>$MYLOG

        cat *_last.maf > read2.maf
        rm *_last.maf *-split.fq

    # Combine Reads
        echo "START Combining Paired End  MAFs [`date`]" >>$MYLOG
            last-pair-probs.py -r read1.maf read2.maf > ${NAME}_last.maf 2>>$MYLOG
        echo "FINISH Combining Paired End  MAFs [`date`]" >>$MYLOG

    # Convert MAF to SAM
        echo "START Converting from MAF to SAM [`date`]" >>$MYLOG
            maf-convert.py sam ${NAME}_last.maf > ${NAME}_last.sam 2>>$MYLOG 
        echo "FINISH Converting to SAM [`date`]" >>$MYLOG

    # PARSE LAST FILES
        echo "START Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG
        perl /scratch/lfs/mcintyre/sem_gof/scripts/parse_last_sam_v2.pl \
            ${NAME}_last.sam \
            ${NAME}_ambig_last.sam \
            ${NAME}_uniq_last.sam \
            > $ALNLOGS/${NAME}_LAST.log

        rm ${NAME}_last.sam 

        echo "FINISH Parsing LAST SAM file into unique and ambiguous [`date`]" >>$MYLOG

}
