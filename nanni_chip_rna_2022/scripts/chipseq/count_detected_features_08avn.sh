#!/bin/sh


## Count and compare detected features between samples

module purge
module load python/3.6

echo "Counting detected features... $(date)"

## Set log file
LOG=${DABG}/count_DAI_${SLURM_JOBID}.log
#echo "DAI LOG FILE: '${LOG}'"
#echo "LOG FOR DAI COUNTING... $(date)" > ${LOG}

## Merge orthologs and GO annotations
if [[ ! -e /ufrc/mcintyre/share/references/dsim_fb202/gene2go_with_terms_from_mel_ortho.csv ]]; then
    python ${SCRIPTS}/merge_go_orthologs_02avn.py \
        -g /ufrc/mcintyre/share/references/dmel_fb617/gene2go_with_terms_02jrbn.csv \
        -o /ufrc/mcintyre/share/etoh_srna/ortholog_files/dmel_orthologs_dsim_fb_2017_04.csv \
        --output /ufrc/mcintyre/share/references/dsim_fb202/gene2go_with_terms_from_mel_ortho.csv \
        >> ${LOG}
fi

## Generate species specific design files
for SPECIES in ${SPECIES1} ${SPECIES2}; do

    ## Make gene list directory
    mkdir -p ${DABG}/${SPECIES}_DAI_gene_lists

    echo "Counting ${SPECIES}..." >> ${LOG}

    ## Get correct fragment fusion and intron files for species
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        FRAGMENT=${FRAGMENT1}
        INTRON=${INTRON1}
        FUSION=${FUSION1}
        GTF=${GTF1}
    else
        FRAGMENT=${FRAGMENT2}
        INTRON=${INTRON2}
        FUSION=${FUSION2}
        GTF=${GTF2}
    fi

    ## Generate temporary fragment fusion and intron files in feature directory
    ## Used in following python script
    cp $(dirname ${FRAGMENT})/$(basename ${FRAGMENT} s_coverage.bed)_annotations.csv \
        ${FEATURES}/${SPECIES}_fragment_xcrpt_annotation.csv
    cp $(dirname ${FUSION})/$(basename ${FUSION} s_coverage.bed)_annotations.csv \
        ${FEATURES}/${SPECIES}_fusion_xcrpt_annotation.csv
    cp $(dirname ${INTRON})/$(basename ${INTRON} .bed).csv \
        ${FEATURES}/${SPECIES}_intron_xcrpt_annotation.csv

    ## Count detected features and generate gene lists
    python ${SCRIPTS}/count_detected_features_08avn.py \
        -i ${DABG}/${SPECIES}_DAI_summary_flags.csv \
        -f ${FEATURES}/${SPECIES}_all_featureID_2_featureType.csv \
        -s ${SPECIES} \
        -d ${ROZ}/temp_${SPECIES}.sql.db \
        -p ${FEATURES} \
        -o ${DABG}/${SPECIES}_chip_features \
        -l ${DABG}/${SPECIES}_DAI_gene_lists \
        >> ${LOG}
    ## Remove temporary flag files
    rm ${ROZ}/temp_${SPECIES}.sql.db

    ## Calculate Cohen's Kappa agreements
    python ${SCRIPTS}/chipseq_kappa_agreement.py \
        -i ${DABG}/${SPECIES}_chip_features_flag.csv \
        -o ${DABG}/${SPECIES}_chip_features \
        >> ${LOG}

    ## Subset flag file for features of interest
    awk -F "," 'NR==1 || $1=="5UTR" || $1=="3UTR" || $1=="TSS300bpWindow" || \
        $1=="fragment" || $1=="intron"  || $1=="intergenic"' \
        ${DABG}/${SPECIES}_chip_features_flag.csv \
        > ${DABG}/${SPECIES}_chip_5U_3U_TSS_frag_intr_inter_flag.csv
    awk -F "," 'NR==1 || $1=="5UTR" || $1=="3UTR" || $1=="TSS300bpWindow" || \
        $1=="fusion" || $1=="intron"  || $1=="intergenic"' \
        ${DABG}/${SPECIES}_chip_features_flag.csv \
        > ${DABG}/${SPECIES}_chip_5U_3U_TSS_fus_intr_inter_flag.csv

    ## Get list of all FBgn in each featureType
#    for feature in 5UTR 3UTR TSS1kbWindow TSS300bpWindow; do
#        awk -F "," 'NR!=1{print $4}' \
#            ${FEATURES}/${SPECIES}_${feature}_xcrpt_annotation.csv | sort | uniq \
#            > ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_${feature}_all.txt
#    done

    ## Remove temporary fragment and intron files
    rm ${FEATURES}/${SPECIES}_fragment_xcrpt_annotation.csv
    rm ${FEATURES}/${SPECIES}_fusion_xcrpt_annotation.csv
    rm ${FEATURES}/${SPECIES}_intron_xcrpt_annotation.csv

    ## Add gene symbol information to TSS300bpWindow lists
    ## (for figure prep)
#    for file in $(ls ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_*_only.txt); do
#        OUTNAME=$(dirname ${file})/$(basename ${file} .txt).geneSymbol_strand.csv
#        python ${SCRIPTS}/get_geneSymbol_strand_chrom_02avn.py \
#            -i ${file} \
#            -g ${GTF} \
#            -c \
#            -o ${OUTNAME} \
#            > $(dirname ${file})/$(basename ${file} .txt).x_vs_auto_counts.txt
#    done
#    echo "Gene counts for all feature types:
#$(wc -l ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_*_all.txt)
#" >> ${LOG}

    ## Make file for JMP genomics with flags for each list
    ## Add gene_id header to main file of all TSS300bpWindow
#    awk '{if(NR==1){print "gene_id\n"$0}else{print $0}}' \
#        ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_all.txt \
#        > ${ROZ}/temp_TSS300bpWindow_all.txt
#    MAIN=${ROZ}/temp_TSS300bpWindow_all.txt
#    echo "${SPECIES} gene list log..." > ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.log
#    for lst in f_both_antibody K27_both_sexes K27_f_only K27_m_only K4_both_sexes K4_f_only K4_m_only m_both_antibody; do
#        echo ${lst} >> ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.log
#        ## Add gene_id to main and list
#        awk '{if(NR==1){print "gene_id\n"$0}else{print $0}}' \
#            ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_${lst}.txt \
#            > ${ROZ}/temp_TSS300bpWindow_${lst}.txt
#        python ${SCRIPTS}/add_flag_gene.py \
#            -m ${MAIN} \
#            -i ${ROZ}/temp_TSS300bpWindow_${lst}.txt \
#            -f ${lst} \
#            -o ${ROZ}/temp_TSS300bpWindow_add_${lst}.txt \
#            >> ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.log
#        MAIN=${ROZ}/temp_TSS300bpWindow_add_${lst}.txt
#    done
#    mv ${MAIN} ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.txt
#    rm ${ROZ}/temp_TSS300bpWindow_*

    ## Run full set of GO annotations (including those without a sim ortholog)
#    if [[ ${SPECIES} == "mel" ]]; then
#        echo "Full set of GO annotations in mel..." >> ${LOG}
#        python ${SCRIPTS}/merge_GO_ID.py \
#            -f ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.txt \
#            -g /ufrc/mcintyre/share/references/dmel_fb617/gene2go_with_terms_02jrbn.csv \
#            --name-flag gene_id \
#            --name-go FBgn \
#            -o ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags_allGO \
#            >> ${LOG}
#    fi
    ## Run single unique pairs of orthologs with GO IDs (no paralogs)
#    echo "${SPECIES} from only single unique pairs of orthologous genes with GO annotations..." >> ${LOG}
#    python ${SCRIPTS}/merge_GO_ID.py \
#        -f ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags.txt \
#        -g /ufrc/mcintyre/share/references/dsim_fb202/gene2go_with_terms_from_mel_ortho.csv \
#        --name-flag gene_id \
#        --name-go ${SPECIES}_geneID \
#        -o ${DABG}/${SPECIES}_DAI_gene_lists/${SPECIES}_TSS300bpWindow_flags_orthoGO \
#        >> ${LOG}
done
