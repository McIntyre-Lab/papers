#!/bin/bash

## Get info from design file
PROJ=/blue/mcintyre/share/maize_ainsworth
DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv

## Make count table
TABLE=$PROJ/MaizeWang2018_transcriptome_eval/MaizeWang2018_consolidated_transcriptome_counts.csv
echo "sample,num_transcript,num_gene" > ${TABLE}

for LINE in $(head -6 ${DESIGN_FILE}); do

    SPECIES=$(echo ${LINE} | awk -F "," '{print $2}')
    TISSUE=$(echo ${LINE} | awk -F "," '{print $3}')

    ## Get sampleID (species_tissue)
    SAMPLEID=${SPECIES}_${TISSUE}

    ## Get count files
    CONSOL=$PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome/${SAMPLEID}_consol_transcript_gene_counts.txt
    NOEA=$PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome/${SAMPLEID}_diff_input_ea_output_transcript.csv

    echo "${SAMPLEID}: $(cat ${NOEA}|wc -l) single transcript genes on Mt/Pt"

    awk -v sample=${SAMPLEID} -v mtpt=$(cat ${NOEA}|wc -l) \
        '{if($0~"transcripts represented"){xcrpt=$1+mtpt;}else if($0~"genes represented")\
        {gene=$1+mtpt}}END{print sample","xcrpt","gene}' ${CONSOL} \
        >> ${TABLE}

done
