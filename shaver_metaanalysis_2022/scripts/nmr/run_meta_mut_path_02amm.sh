#!/bash/bin


PROJ=~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022
SCRIPTS=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/scripts/NMR/meta_analysis

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH


mkdir -p $PROJ/nmr/meta_plotting

#python $SCRIPTS/meta_mut_path_loop.py

Rscript $SCRIPTS/heatmap_meta_effect_loop.R

