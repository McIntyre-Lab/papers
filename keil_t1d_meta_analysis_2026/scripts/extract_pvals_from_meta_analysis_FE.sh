#!/bin/bash

##Extract pvalues from meta-analysis summary
PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

##Set paths - overall metanalysis
IND_ALL=${PROJ}/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary
SCRIPTS=${PROJ}/scripts

ROZ_ALL=${IND_ALL}/roz
   mkdir -p ${ROZ_ALL}

##cp files for extraction to roz
for file in ${IND_ALL}/summary_sexMod_*; do
    cp "$file" ${ROZ_ALL}
done

#cp ${IND_ALL}/summary_sexMod_* ${ROZ_ALL}


#Extract pvalues and output CSV - overall meta
python ${SCRIPTS}/extract_pvals_from_meta_analysis_summary.py \
    -d ${ROZ_ALL}\
    -s FE.txt\
    -x 0.05 \
    -o ${IND_ALL}/meta_analysis_pvals_moderator_heterogenity.csv

rm -r ${ROZ_ALL}

#Exons

##Set paths - metanalysis exons only
IND_EXON=${PROJ}/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_exonsOnly

ROZ_EXON=${IND_EXON}/roz
   mkdir -p ${ROZ_EXON}

##cp files for extraction to roz
for file in ${IND_EXON}/summary_sexMod_*; do
    cp "$file" ${ROZ_EXON}
done


#Extract pvalues and output CSV - metanalysis exons only
python ${SCRIPTS}/extract_pvals_from_meta_analysis_summary.py \
    -d ${ROZ_EXON}\
    -s FE.txt\
    -x 0.05 \
    -o ${IND_EXON}/meta_analysis_pvals_moderator_heterogenity_exonsOnly.csv

rm -r ${ROZ_EXON}

#Introns

##Set paths - metanalysis introns only
IND_INTRON=${PROJ}/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_intronsOnly

ROZ_INTRON=${IND_INTRON}/roz
   mkdir -p ${ROZ_INTRON}

##cp files for extraction to roz
for file in ${IND_INTRON}/summary_sexMod_*; do
    cp "$file" ${ROZ_INTRON}
done


#Extract pvalues and output CSV - metanalysis introns only
python ${SCRIPTS}/extract_pvals_from_meta_analysis_summary.py \
    -d ${ROZ_INTRON}\
    -s FE.txt\
    -x 0.05 \
    -o ${IND_INTRON}/meta_analysis_pvals_moderator_heterogenity_intronsOnly.csv

rm -r ${ROZ_INTRON}
