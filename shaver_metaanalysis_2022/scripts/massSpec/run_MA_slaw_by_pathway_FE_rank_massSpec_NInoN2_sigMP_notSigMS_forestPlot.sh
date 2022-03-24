#!/bash/bin

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH


PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPT=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

READY=$PROJ/SLAW_UGA_Output/analysisReady_datasets

DESIGN=$PROJ/design_files/dsgn_GT_${VARDIR}_pairs_slaw.tsv
    #sampleID        sample  strain  batch   sampleType      newBatch        sampleType_batch        strain_w_poolPD pair    set
    #b1_aos_122      aos_122 RB2011  1       strain  b1      strain_b1       RB2011  1       set1
    #b1_aos_127      aos_127 CX11314 1       strain  b1      strain_b1       CX11314 1       set1

OUTD=$PROJ/SLAW_UGA_Output/meta_analysis_${VAR}
    mkdir -p $OUTD

ROZ=$PROJ/SLAW_UGA_Output/roz_MA/${VAR}
    mkdir -p $ROZ 


## subset data to only features sig in meta path and not in individual NI strains (omitted N2)
head -n1 $READY/${VAR}_analysisReady_rank_sbys.tsv > $ROZ/plots_${VAR}_rank_sbys.tsv

# convert list to grep format and pull rows from starting dataset
#LIST=`(sed '1d' $PROJ/SLAW_UGA_Output/meta_analysis_${VAR}/list_features_sig_metaPath_NInon2_notSig_metaStrain_${VAR}.tsv | tr '\n' '|' | sed 's/.$//')`
#LIST=$PROJ/SLAW_UGA_Output/meta_analysis_${VAR}/list_features_sig_metaPath_NInon2_notSig_metaStrain_${VAR}.tsv
#grep -E $LIST $READY/${VAR}_analysisReady_rank_sbys.tsv >> $ROZ/plots_${VAR}_rank_sbys.tsv

cat $PROJ/SLAW_UGA_Output/meta_analysis_${VAR}/list_features_sig_metaPath_NInon2_notSig_metaStrain_${VAR}.tsv | while read line
do
    cat $READY/${VAR}_analysisReady_rank_sbys.tsv | grep $line$'\t' >> $ROZ/plots_${VAR}_rank_sbys.tsv
done

## run MA
PATHWAY="NInoN2"
awk '$3=="CB4856" || $3=="CX11314" || $3=="DL238" || $3=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
    --wide $ROZ/plots_${VAR}_rank_sbys.tsv \
    --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
    -id uniqueID \
    --study "batch" \
    --treatment "strain" \
    --contrast "CB4856,CX11314,DL238,PD1074" \
    --report $OUTD/MA_FE_path_NInoN2_sigMP_notSigMS_${VAR}_report.tsv \
    -o $OUTD/MA_FE_path_NInoN2_sigMP_notSigMS_${VAR}.tsv \
    --forest $OUTD/MA_FE_path_NInoN2_sigMP_notSigMS_${VAR} \
    --model FE \
    --effectSize SMD \
    --cmMethod UB
