#!/bin/bash

###############################
### Collect unique features ###
###############################


module purge
module load bedtools/2.28.0
module load samtools/1.9
module load python/3.6

## Set output directory
export FEATuniq=$PROJ/dros_annotation_pipeline
    mkdir -p ${FEATuniq}

## Run loop for each species
for SPECIES in ${SPECIES1} ${SPECIES2}; do
    echo "${SPECIES}...
"
    ## Set temporary directory
    export ROZ=${FEATuniq}/temp_${SPECIES}
        mkdir -p ${ROZ}

    ## Set reference FASTA, GTF, and EA annotations
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        FASTA=${FASTA1}
        GTF=${GTF1}
        FRAGMENT=${EA1}_exon_fragments_coverage.bed
        FRAGMENTannot=${EA1}_exon_fragment_annotations.csv
        FUSION=${EA1}_fusions_coverage.bed
        FUSIONannot=${EA1}_fusion_annotations.csv
        INTRON=${EA1}_introns_from_fusions.bed
        INTRONannot=${EA1}_introns_from_fusions.csv
    else
        FASTA=${FASTA2}
       	GTF=${GTF2}
       	FRAGMENT=${EA2}_exon_fragments_coverage.bed
        FRAGMENTannot=${EA2}_exon_fragment_annotations.csv
        FUSION=${EA2}_fusions_coverage.bed
        FUSIONannot=${EA2}_fusion_annotations.csv
       	INTRON=${EA2}_introns_from_fusions.bed
        INTRONannot=${EA2}_introns_from_fusions.csv
    fi

    ## Generate chrom.sizes file for each reference
    echo -e "\tGenerate chrom.sizes..."
    source ${SCRIPTS}/make_chrom_sizes.sh

    ## Obtain feature coordinates from reference annotation
    ## for all annotated 5'UTR, 3'UTR, and TSS (1kb region centered on start of transcript)
    echo -e "\tGet 5UTR, 3UTR, and TSS1kbWindow coordinates..."
    python ${SCRIPTS}/get_features_06avn.py \
        -i ${GTF} \
        -g ${CHROM} \
        -o ${ROZ}

    ## Get unique values for 5'UTR, 3'UTR, and TSS1kbWindow
    echo -e "\tGet unique values for 5UTR, 3UTR, and TSS1kbWindow..."
    for feature in 5UTR 3UTR TSS1kbWindow; do
        python ${SCRIPTS}/unique_features_05avn.py \
            -b ${ROZ}/all_${feature}.bed \
            -o ${FEATuniq}/${SPECIES}_${feature}
    done

    ## Get set of unique TSS300bpWindow
    echo -e "\tGet unique coordinates for TSS300bpWindow..."
    python ${SCRIPTS}/make_unique_TSS300bpWindow_03avn.py \
        -i ${GTF} \
        -g ${CHROM} \
        -o ${FEATuniq}/${SPECIES}_TSS300bpWindow


    ## Identify intergenic regions
    echo -e "\tIdentify intergenic regions..."
    source ${SCRIPTS}/intergenic_03avn.sh


    ## Combine and summarize unique feature counts
    echo -e "\tCombine and summarize unique features..."
    source ${SCRIPTS}/combine_summ_unique_feat_02avn.sh


    ## Connect featureID to featureType
    echo -e "\tGenerate featureID_2_featureType..."
    source ${SCRIPTS}/make_feature_type_file_02avn.sh


    ## Connect featureID to all associated FBgn
    echo -e "\tGenerate featureID_2_fbgn..."
    source ${SCRIPTS}/make_featureID_2_fbgn_02avn.sh


    ## Remove temporary directory
    #rm -r ${ROZ}

done

## Count regions
## Number of genes, transcripts, fragments, fusions,
##     introns, and intergenic regions in each species
echo -e "
Count regions in both species..."
source ${SCRIPTS}/count_genes_transcripts_02avn.sh

