#!/bash/bin


PROJ=~/mclab/SHARE/McIntyre_Lab/CID/sweet16

SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH

for i in cdcl3 d2o
do

DATA=$PROJ/nmr/meta_analysis

python $SCRIPTS/upset_plot.py \
    --wide $DATA/MA_byStrain_sig_4_upset_${i}.csv \
    --design $DATA/upset_dsgn_${i}.csv \
    --uniqID featureID \
    --title upsetPlot_byStrain_sig_${i}_min3 \
    --outD $DATA/
done
