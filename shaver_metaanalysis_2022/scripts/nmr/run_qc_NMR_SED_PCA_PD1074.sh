#!/bin/bash


#bash
#source activate /TB14/TB14/gait_gm/galaxy/database/dependencies/_conda/envs/__gaitGm@1.0.1
export PATH=/TB14/TB14/conda_envs/gait-gm/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts


## SED and PCA on NMR sweet16 
    ## only PD1074 in more that 1 batch


INPUT=$PROJ/data_from_group/nmr_sweet16
QC=$INPUT/qc_PD1074
    mkdir -p $QC

for i in CDCL3 D2O
do
for j in all mult
do

## SED
    echo "SED ${i}
    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/ph_NMR_${i}_shortNames.tsv \
        --design $INPUT/dsgn_NMR_${i}_${j}_PD1074.txt \
        --ID ppm \
        --figure $QC/SED_figure_NMR_${i}_${j}_PD1074.pdf \
        --SEDtoMean $QC/SEDtoMean_NMR_${i}_${j}_PD1074.tsv \
        --SEDpairwise $QC/SEDpairwise_NMR_${i}_${j}_PD1074.tsv \
        --per 0.8


## PCA
    echo "PCA $i
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/ph_NMR_${i}_shortNames.tsv \
        --design $INPUT/dsgn_NMR_${i}_${j}_PD1074.txt \
        --ID ppm \
        --load_out $QC/PCA_loading_NMR_${i}_${j}_PD1074.tsv \
        --score_out $QC/PCA_scores_NMR_${i}_${j}_PD1074.tsv \
        --summary_out $QC/PCA_summary_NMR_${i}_${j}_PD1074.tsv \
        --figure $QC/PCA_figure_pos_NMR_${i}_${j}_PD1074.pdf
done
done
