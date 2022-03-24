#!/bash/bin


PROJ=~/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=$PROJ/scripts/massSpec/meta_analysis

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH


mkdir -p $PROJ/SLAW_UGA_Output/comparing_models


#python $SCRIPTS/meta_mut_path_loop.py

Rscript $SCRIPTS/heatmap_meta_effect_loop.R

