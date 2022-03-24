#!/bin/bash


#bash
#source activate /TB14/TB14/gait_gm/galaxy/database/dependencies/_conda/envs/__gaitGm@1.0.1
export PATH=/TB14/TB14/conda_envs/gait-gm/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts


## PLS-DA NMR sweet16 
    
INPUT=$PROJ/data_from_group/nmr_sweet16
QC=$INPUT/qc
    mkdir -p $QC

## compare PD1074 to all others
for i in CDCL3 D2O
do

for j in AUM2073 CB4856 CX11314 DL238 KJ550 N2 RB2011 RB2055 RB2550 UGT49 UGT60 VC1265 VC2524
do

## PLS-DA
echo "PLS-DA ${i}
"

python $SCRIPTS/partial_least_squares.py \
    --input $INPUT/ph_NMR_${i}_shortNames.tsv \
    --design $INPUT/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --toCompare PD1074,${j} \
    -cv none \
    --nComp 2 \
    --outScores $QC/PLSDA_loadings_NMR_${i}_PD1074_vs_${j}.tsv \
    --outWeights $QC/PLSDA_weights_NMR_${i}_PD1074_vs_${j}.tsv \
    --outClassification $QC/PLSDA_class_NMR_${i}_PD1074_vs_${j}.tsv \
    --outClassificationAccuracy $QC/PLSDA_class_accuracy_NMR_${i}_PD1074_vs_${j}.tsv \
    --figure $QC/PLSDA_figure_NMR_${i}_PD1074_vs_${j}.pdf

done
done 

