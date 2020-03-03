#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rsemStar
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/rsemStar_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
###SBATCH --array=1-120
#SBATCH --array=1

## Using looping of the array task runs
## df (no header) has 360 rows
## want 5 rows per array task
        ## split into X array tasks:  120/X = 5 rows per array task, X=24


#Set the number of runs that each SLURM task should do
##PER_TASK=5
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    echo "This is SLURM task $SLURM_ARRAY_TASK_ID, run number $RUN"

    module purge
    module load R
    module load bowtie/0.12.9
    module load rsem/1.2.28

    ### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    REF=$PROJ/pacbio_analysis/references

    RSEMREF=$REF/genome_rsem_index_STAR

     DATAOUT=$PROJ/rsemStar_dataPrep_junctions
        mkdir -p $DATAOUT

    ## DF: maize_rnaseq_RG_sampleCode_list.csv
        #      1 sampleID,RAPiD_Genomics_Sample_Code
        #      2 B73_P1_C1_Amb,UFE_032201_P003_WE01
        #      3 B73_P1_C2_Amb,UFE_032201_P003_WA08
        #      4 B73_P1_C4_Ele,UFE_032201_P003_WE06
        #      5 B73_P1_C5_Ele,UFE_032201_P003_WA02
        #      6 B73_P1_C6_Amb,UFE_032201_P003_WD09
        #      7 B73_P1_C7_Ele,UFE_032201_P003_WA05
        #      8 B73_P2_C1_Amb,UFE_032201_P003_WA10
        #      9 B73_P2_C2_Amb,UFE_032201_P003_WD10

    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_rnaseq_RG_sampleCode_list.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    RG=${ARRAY[1]}

    $ROZ=$PROJ/tmp_rozStar/${SAMPLE}
        mkdir -p $ROZ

    ## cat together lanes and reads --> aligning raw reads SE
    cat $DATAIN/RAPiD-Genomics_F086_${RG}_* > $ROZ/$SAMPLE}_all.fq
    READS=$ROZ/$SAMPLE}_all.fq


    ## (1) SE alignment of catted raw reads to genome ref with star

    ## default aligner is bowtie1, or request --bowtie2  --star
    ## note jeremy used following options:
	##        --fragment-length-mean 80 \
	##        --fragment-length-sd 50 \

    OUTSE=$OUTDIR/${SAMPLE}_SE
        mkdir -p $OUTSE

: <<'END'
    cat $INPUT/${SAMPLE}_*_bbmerge_min*_distinct.fq > $ROZ/${SAMPLE}_bbmerge.fq
    READS=$ROZ/${SAMPLE}_bbmerge.fq

    rsem-calculate-expression -p 4 \
        --temporary-folder ${ROZ} \
        --calc-ci \
        --estimate-rspd \
        ${READS} ${RSEMREF} ${OUTSE}

END
    ## (2) PE alignment of unmerged reads to reference PB transcriptome ==>  {SAMPLE}_techRep${TR}_R1/2_unmerge_min${MIN}_distinct.fq

    OUTPE=$OUTDIR/${SAMPLE}_PE
        mkdir -p $OUTPE

: <<'END'
    cat $INPUT/${SAMPLE}_*_R1_unmerge_min*_distinct.fq > $ROZ/${SAMPLE}_R1_unmerge.fq
    cat $INPUT/${SAMPLE}_*_R2_unmerge_min*_distinct.fq > $ROZ/${SAMPLE}_R2_unmerge.fq
    R1=$ROZ/${SAMPLE}_R1_unmerge.fq
    R2=$ROZ/${SAMPLE}_R2_unmerge.fq

    rsem-calculate-expression -p 4 \
        --paired-end \
        --temporary-folder ${ROZ} \
        --calc-ci \
        --estimate-rspd \
        $R1 $R2 \
        ${READS} ${RSEMREF} ${OUTPE}
END


    ### sum SE and PE matrices (cat together and use groupby(sum) transcript

    ### prep matrix - gene and transcript
#    SECNT=$(ls $OUTDIR/*_SE.isoforms.results | wc -l) 
#    PECNT=$(ls $OUTDIR/*_PE.isoforms.results | wc -l)    

    for i in SE
    do
        cd $OUTDIR
        for samp in *_${i}.isoforms.results
        do
          echo ${samp}
          cut -f1,5 ${samp} > $ROZ/${samp}".expected_count"
          s=$(echo "${samp}" | sed 's/_'"${i}"'.isoforms.results//')
          echo ${s}
          sed s/'expected_count'/${s}/g $ROZ/${samp}".expected_count" > $ROZ/${s}_${i}.exp
        done

        ## paste all output together and keep transcriptID and expected count for each sample
        paste $ROZ/*_${i}.exp > $ROZ/${i}_all.exp

        COLS=$(awk '{ cols = "1"; \
            for (x=1; x<=NF; x++) \
            {if (x % 2 ==0) \
                { cols = cols","x; \
                } } } END {print cols}' $ROZ/${i}_all.exp) 
        echo $COLS

        cut -f $COLS $ROZ/${i}_all.exp > $OUTDIR/${i}_expression_matrix.pre-filter.expected_count.txt
   done

##    paste $ROZ/*_${i}'.exp' | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48 > expression_matrix.pre-filter.exp_count.txt


##    paste $OUT/${SAMPLE}/rsem.genes.results | tail -n+2 | cut -f1,5,12,19,26 > $OUT/${SAMPLE}_genes.rsem.txt
#    paste $OUT/${SAMPLE}/rsem.isoforms.results | tail -n+2	| cut -f1,5,12,19,26 > $OUT/${SAMPLE}_isoforms.rsem.txt

##rm -r $ROZ

done
