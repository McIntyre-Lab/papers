#!/bin/bash

# for samples NOT b73, identify transcripts that do NOT map to ref B73 but DO map to ref Mo17


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
END_NUM=9

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_noB73_samples.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    GENO=${ARRAY[1]}
    TRT=${ARRAY[2]}
    ID=${ARRAY[3]}

    SAMPLE=${ID}_${GENO}_${TRT}

    LIST=$INPUT/list_${SAMPLE}.ignored_unmapped_ids.b73_ref.txt
    SAM=/home/ammorse/TB14/maize_ainsworth/mapping_minimap_mo17_yan/${SAMPLE}.polished.all.hq.mapped.sorted.sam
    INDEX=/home/ammorse/TB14/maize_ainsworth/mapping_minimap_mo17_yan/Zm-Mo17-REFERENCE-YAN-1.0.fasta.fai


    echo -e "num_NOT_ALN_2_Mo17_yan	sampleID" > $ROZ/no_header.txt
    echo -e "num_ALN_2_Mo17_yan	sampleID" > $ROZ/yes_header.txt

    while read TRANSCRIPTS
    do
        samtools view -f 4 -t $INDEX $SAM | grep "${TRANSCRIPTS}" | wc -l >> $ROZ/cnts_no_${SAMPLE}.txt
        samtools view -F 4 -t $INDEX $SAM | grep "${TRANSCRIPTS}" | wc -l >> $ROZ/cnts_yes_${SAMPLE}.txt
    done < $LIST

    awk -v sample="$SAMPLE" '{sum += $1} END {print sum "\t" sample}' $ROZ/cnts_no_${SAMPLE}.txt >> $ROZ/${SAMPLE}_sum_no.txt
    awk -v sample="$SAMPLE" '{sum += $1} END {print sum "\t" sample}' $ROZ/cnts_yes_${SAMPLE}.txt >> $ROZ/${SAMPLE}_sum_yes.txt

done

cat $ROZ/no_header.txt $ROZ/*_sum_no.txt > $OUTPUT/no_2_Mo17_yan.txt
cat $ROZ/yes_header.txt $ROZ/*_sum_yes.txt > $OUTPUT/yes_2_Mo17_yan.txt

rm -r $ROZ
