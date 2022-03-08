#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=extract
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/extract_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --array=1-11

module load python/3.8
module load seqtk

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth

for i in b73 mo17_yan mo17_cau
do

    OUT=$PROJ/compare_b73_2_mo17/uncollapsed_${i}_IsoSeq_reads
        mkdir -p $OUT

    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    GENO=${ARRAY[1]}
    TRT=${ARRAY[2]}    
    ID=${ARRAY[3]}

    SAMPLE=${ID}_${GENO}_${TRT}

    INPUT=$PROJ/compare_b73_2_mo17/tofu2_${i}/${SAMPLE}

    ## DONE PREVIOUSLY
    ## counts ==> all, unmapped, low coverage, low identity
#    wc -l $INPUT/${SAMPLE}.ignored_ids.txt | awk -v sample=$SAMPLE '{print "all" "\t" sample "\t" $1}' > $OUT/cnt_ignored_all_${SAMPLE}_ref_${i}.txt
#    grep "Unmapped" $INPUT/${SAMPLE}.ignored_ids.txt | wc -l | awk -v sample=$SAMPLE '{print "unmapped" "\t" sample "\t" $0}' > $OUT/cnt_ignored_unmapped_${SAMPLE}_ref_${i}.txt
#    grep "Coverage" $INPUT/${SAMPLE}.ignored_ids.txt | wc -l | awk -v sample=$SAMPLE '{print "low_coverage" "\t" sample "\t" $0}' > $OUT/cnt_ignored_coverLow_${SAMPLE}_ref_${i}.txt    
#    grep "Identity" $INPUT/${SAMPLE}.ignored_ids.txt | wc -l | awk -v sample=$SAMPLE '{print "low_identity" "\t" sample "\t" $0}' > $OUT/cnt_ignored_identityLow_${SAMPLE}_ref_${i}.txt
    ## create list of 'ignored', single column only
#    awk '{print $1}' $INPUT/${SAMPLE}.ignored_ids.txt > $OUT/list_${SAMPLE}.ignored_ids.${i}_ref.txt

    ## only low identity reads
    ## create list of 'ignored and low identity', single column only
    grep "Identity" $INPUT/${SAMPLE}.ignored_ids.txt | awk '{print $1}' > $OUT/list_${SAMPLE}.ignored_low_ident_ids.${i}_ref.txt
    
    LIST=$OUT/list_${SAMPLE}.ignored_low_ident_ids.${i}_ref.txt
    FA=$PROJ/compare_b73_2_mo17/mapping_minimap2_${i}/${GENO}/${TRT}/${SAMPLE}.polished.all.hq.fasta
    FQ=$PROJ/compare_b73_2_mo17/mapping_minimap2_${i}/${GENO}/${TRT}/${SAMPLE}.polished.all.hq.fastq
    echo "
	list of ignored low identity IDs:  $LIST
        input fasta file to search:  $FA
        input fastq file to search: $FQ
        "
    ## Get FA subset
#    python $PROJ/scripts/Common_extracFeatureFromList.py \
#        --list ${LIST} \
#        --type fasta \
#        --input ${FA} \
#        --out $OUT/${SAMPLE}_unCollapsed_low_ident_${i}_ref.fa

    ## Get FQ subset
    seqtk subseq ${FQ} ${LIST} > $OUT/${SAMPLE}_unCollapsed_low_ident_${i}_ref.fq

done
