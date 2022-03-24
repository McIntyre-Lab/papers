#!/bin/sh

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts



##  for features sig in at least 1 strain

for TYPE in cdcl3 d2o
do

## rotate data - metabolites across top
datamash transpose < $PROJ/nmr/meta_analysis/MA_byStrain_sig_4_heatmap_${TYPE}.tsv > $PROJ/nmr/meta_analysis/MA_byStrain_sig_4_heatmap_${TYPE}_rotated.tsv

## design for rotated data
head -n1 $PROJ/nmr/meta_analysis/MA_byStrain_sig_4_heatmap_${TYPE}_rotated.tsv | tr '\t' '\n' | sed 's/featureID/sampleID/' > $PROJ/nmr/meta_analysis/dsgn_effectSizes_4_heatmap_rotated.tsv




#MMC metabolites across  top
python $SCRIPTS/modulated_modularity_clustering.py \
    --input $PROJ/nmr/meta_analysis/MA_byStrain_sig_4_heatmap_${TYPE}_rotated.tsv \
    --design $PROJ/nmr/meta_analysis/dsgn_effectSizes_4_heatmap_rotated.tsv \
    --ID featureID \
    --correlation pearson \
    --sigmaLow 0.05 \
    --sigmaHigh 0.50 \
    --sigmaNum 451 \
    --figure $PROJ/nmr/meta_analysis/MMC_heatmap_MA_byStrain_sig_${TYPE}.pdf \
    --out $PROJ/nmr/meta_analysis/MMC_output_MA_byStrain_sig_${TYPE}.tsv

done
