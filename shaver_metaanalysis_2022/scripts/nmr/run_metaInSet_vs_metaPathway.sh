#!/bash/bin




export PATH=/TB14/TB14/conda_envs/gait-gm/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=$PROJ/scripts/NMR/meta_analysis


mkdir -p $PROJ/nmr/meta_analysis/comparing_models


python $SCRIPTS/metaInSet_vs_metaPathway_02amm.py

Rscript $SCRIPTS/heatmap_meta_effect_02amm.R




