#!/bash/bin


PROJ=~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec

SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts

DATA=$PROJ/upset_plots

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH



for i in rp_neg rp_pos hilic_pos
do
    for j in TCA UGT NI NInoN2 Combo
    do

## sig
python $SCRIPTS/upset_plot.py \
    --wide $DATA/${i}_${j}_4_upset_plot.csv \
    --design $DATA/upset_plot_design_file_${j}_sig.csv \
    --uniqID featureID \
    --title upset_plot_sig_${i}_${j}_sig \
    --outD $DATA/

## direction
python $SCRIPTS/upset_plot.py \
    --wide $DATA/${i}_${j}_4_upset_plot.csv \
    --design $DATA/upset_plot_design_file_${j}_ge_pd1074.csv \
    --uniqID featureID \
    --title upset_plot_sig_${i}_${j}_dir_ge_pd1074 \
    --outD $DATA/

    done
done
