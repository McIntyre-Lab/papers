#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=star
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/star-%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --cpus-per-task=8
#SBATCH --array=2-121


## df has header (2-121)

### star to generate junction alignments
	## readLength is 150bp
	## length on either side of junction is:  (RL/2 + 13)
        ## juncLength is:  (150/2 + 13) *2 = 176

## aligning RAW reads

    module load star

    ### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    NUMPROC=8
    REFDIR=$PROJ/pacbio_analysis/references/star_index_mo17
    DATAIN=$PROJ/original_data/RNAseq/Plate3-4

    DATAOUT=$PROJ/star_dataPrep_junctions_Mo17
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

    ROZ=$PROJ/rozStar/roz_${SAMPLE}
        mkdir -p $ROZ

    ## cat together lanes  --> aligning raw reads PE
    cat $DATAIN/RAPiD-Genomics_F086_${RG}_*_R1_*.fastq.gz > $ROZ/${SAMPLE}_R1.fq.gz
    cat $DATAIN/RAPiD-Genomics_F086_${RG}_*_R2_*.fastq.gz > $ROZ/${SAMPLE}_R2.fq.gz
    READ1=$ROZ/${SAMPLE}_R1.fq.gz
    READ2=$ROZ/${SAMPLE}_R2.fq.gz

#    READ1=$(gunzip -c $ROZ/${SAMPLE}_R1.fq.gz)
#    READ2=$(gunzip -c $ROZ/${SAMPLE}_R2.fq.gz)

echo "$(date)
    aligning raw reads PE
 
    reads are :  $READ1 
   		 $READ2
    
    reduce spurious junctions (--outFilterType BySJout)
    min overhang for unannotated junctions (--alignSJoverhangMin 8)
    min overhang for annotated junctions (--alignSJDBoverhangMin 1)
    
"

    STAR --runThreadN ${NUMPROC} \
        --genomeDir ${REFDIR} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${READ1} ${READ2} \
        --outFilterType BySJout \
        --outFileNamePrefix $DATAOUT/${SAMPLE}_ \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterMismatchNmax 999 \
         --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 1 \
        --genomeLoad NoSharedMemory \
        --outSAMtype BAM Unsorted

echo "done"

rm -r $ROZ
