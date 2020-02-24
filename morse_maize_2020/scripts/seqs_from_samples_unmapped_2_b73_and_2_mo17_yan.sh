#!/bin/bash

# for samples NOT b73 and NOT Mo17,  extract sequences that do NOT map to b73 NOR to mo17 ya


### Set Directories
PROJ=/home/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
INPUT=$PROJ/evidence_genotype_specific_seqs/uncollapsed_IsoSeq_reads
OUTPUT=$PROJ/evidence_genotype_specific_seqs

ROZ=$OUTPUT/roz
    mkdir -p $ROZ

#LOG=$PROJ/evidence_genotype_specific_seqs/log_file.log

## list mo17 transcripts that do NOT map to b73
#	list_19_mo17_amb.ignored_unmapped_ids.b73_ref.txt
#	list_21-2_mo17_oz.ignored_unmapped_ids.b73_ref.txt
#	list_21_mo17_oz.ignored_unmapped_ids.b73_ref.txt

## do transcripts in the above lists map?
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/amb/19_mo17_amb.polished.all.hq.fasta
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/oz/21-2_mo17_oz.polished.all.hq.fasta
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/oz/21_mo17_oz.polished.all.hq.fasta

START_NUM=1
END_NUM=6

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_noB73_noMo17_samples.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    GENO=${ARRAY[1]}
    TRT=${ARRAY[2]}
    ID=${ARRAY[3]}

    SAMPLE=${ID}_${GENO}_${TRT}

    LIST=$INPUT/list_${SAMPLE}.ignored_unmapped_ids.b73_ref.txt
    SAM=/home/ammorse/TB14/maize_ainsworth/mapping_minimap_mo17_yan/${SAMPLE}.polished.all.hq.mapped.sorted.sam
    INDEX=/home/ammorse/TB14/maize_ainsworth/mapping_minimap_mo17_yan/Zm-Mo17-REFERENCE-YAN-1.0.fasta.fai


    ## extract sequences that do NOT map to b73  NOR to mo17 yan
    while read TRANSCRIPTS
    do
        samtools view -f 4 -t $INDEX $SAM | grep "${TRANSCRIPTS}" >> $ROZ/no_${SAMPLE}.txt
    done < $LIST

    awk -v sample="$SAMPLE" '{ printf ">%s\n%s\n", $1, $10 }' $ROZ/no_${SAMPLE}.txt > $OUTPUT/seqs_no_aln_${SAMPLE}.fa

done

rm -r $ROZ
