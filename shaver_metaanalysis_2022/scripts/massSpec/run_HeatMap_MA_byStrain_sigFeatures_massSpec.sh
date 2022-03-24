#!/bin/bash


#bash
#source activate /TB14/TB14/gait_gm/galaxy/database/dependencies/_conda/envs/__gaitGm@1.0.1
export PATH=/TB14/TB14/conda_envs/gait-gm/bin:$PATH


PROJ=~/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts


### heatmap effect sizes from meta-strain all but N2
##  for features sig in at least 1 strain

for TYPE in rp_neg rp_pos hilic_pos
do

    python $SCRIPTS/hierarchical_clustering_heatmap.py \
        --input $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MA_byStrain_sig_4_heatmap.tsv \
        --design $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/dsgn_effectSizes_4_heatmap.tsv \
        --uniqID featureID \
        --dendogram \
        --fig $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/HeatMap_MA_byStrain_sig_${TYPE}.pdf

done

