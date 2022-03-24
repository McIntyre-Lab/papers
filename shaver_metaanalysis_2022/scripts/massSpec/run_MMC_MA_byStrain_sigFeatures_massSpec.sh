#!/bin/sh

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts


##  for features sig in at least 1 strain

for TYPE in rp_neg rp_pos hilic_pos
do

## rotate data - metabolites across top
datamash transpose < $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MA_byStrain_sig_4_heatmap.tsv > $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MA_byStrain_sig_4_heatmap_rotated.tsv

## design for rotated data
head -n1 $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MA_byStrain_sig_4_heatmap_rotated.tsv | tr '\t' '\n' | sed 's/featureID/sampleID/' > $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/dsgn_effectSizes_4_heatmap_rotated.tsv

#MMC effect size across top
python $SCRIPTS/modulated_modularity_clustering.py \
    --input $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MA_byStrain_sig_4_heatmap_rotated.tsv \
    --design $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/dsgn_effectSizes_4_heatmap_rotated.tsv \
    --ID featureID \
    --correlation pearson \
    --sigmaLow 0.05 \
    --sigmaHigh 0.50 \
    --sigmaNum 451 \
    --figure $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MMC_heatmap_MA_byStrain_sig_${TYPE}.pdf \
    --out $PROJ/SLAW_UGA_Output/meta_analysis_${TYPE}/MMC_output_MA_byStrain_sig_${TYPE}.tsv

done

