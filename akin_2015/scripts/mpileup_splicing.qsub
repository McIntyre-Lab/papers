#!/bin/bash
#PBS -M jrbnewman@ufl.edu
#PBS -m n
#PBS -r n
#PBS -q bio
#PBS -j oe
#PBS -o /scratch/lfs/sugrue/scripts/PBS_LOGS/mpileup_splicing_bt
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=4gb
#PBS -t 1-64

module load samtools/0.1.19

#### Set Directories
    PROJ=/scratch/lfs/sugrue

#### Because I am using an PBS Array I am pulling LINE:MV:REP from an external CSV with all possible combinations
    DESIGN_FILE=$PROJ/design_files/sugrue_file_list.csv
    DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    C=${ARRAY[0]}
    NUM=${ARRAY[1]}
    SAMP=${ARRAY[2]}
    LANE=${ARRAY[3]}
    READ=${ARRAY[4]}
    BIN=${ARRAY[5]}

    NAME=${C}-${NUM}_${SAMP}_${LANE}_${READ}_${BIN}

    # Create OUTPUT directory if needed.
        #OUTPUT=$PROJ/mpileup_splicing_bt
        OUTPUT=$PROJ/mpileup_splicing
        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

        BAMOUT=$PROJ/aln_hg19_bam_splicing
        if [ ! -e $BAMOUT ]; then mkdir -p $BAMOUT; fi

    # Create LOG directory and start log
        LOGS=$OUTPUT/logs 
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        MYLOG=$LOGS/${NAME}.log
        printf "`date` $NAME PBS_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}"

#### References
    REF=/scratch/lfs/sugrue/references/splicing/hg19_splicing_catalogue_176bp.fa

    #if [ ! -e $REF.fai ]
    #then
    #    samtools faidx $REF
    #fi

#### Convert SAM to BAM and make mpileups
    #SAM=$PROJ/bwa_mem_aln_splice/${NAME}.sam 
    SAM=$PROJ/aln_splicing_se/${NAME}.sam
    BAM=$BAMOUT/${NAME}

    printf "<-------------------- Convert SAM to BAM -------------------->" >> "${MYLOG}"
    echo `date`": Starting SAM to BAM conversion" >> "${MYLOG}"
    samtools view -ut $REF.fai -o $BAM.bam $SAM 2>> "${MYLOG}" 
    samtools sort -m 1000000000 $BAM.bam $BAM.sorted 2>> "${MYLOG}"
    samtools index $BAM.sorted.bam >> "${MYLOG}"
    echo `date`": Finished SAM to BAM conversion" >> "${MYLOG}"

#### Make mpielup
    PILEUP=$OUTPUT/${NAME}.mpileup

    printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> "${MYLOG}"
    echo `date`": Generating pileup" >> "${MYLOG}"
    samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> "${MYLOG}"


echo `date`": Finished Script complete" >> "${MYLOG}"
