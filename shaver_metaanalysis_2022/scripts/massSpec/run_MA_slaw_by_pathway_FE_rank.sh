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

PATHWAY="NInoN2"
awk '$3=="CB4856" || $3=="CX11314" || $3=="DL238" || $3=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
    --wide $READY/${VAR}_analysisReady_rank_sbys.tsv  \
    --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
    -id uniqueID \
    --study "batch" \
    --treatment "strain" \
    --contrast "CB4856,CX11314,DL238,PD1074" \
    --report $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_report.tsv \
    -o $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_SMD_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

: <<'END'

## pull out following strain for meta pathway (VC1265=TCA,DL2238=NI,UGT60(VC2512)=UGT
PATHWAY="COMBO"
awk '$3=="VC1265" || $3=="DL238" || $3=="UGT60" || $3=="PD1074" || NR==1' ${DESIGN} > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
    --wide $READY/${VAR}_analysisReady_rank_sbys.tsv \
    --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
    -id uniqueID \
    --study "batch" \
    --treatment "strain" \
    --contrast "VC1265,DL238,UGT60,PD1074" \
    --report $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_report.tsv \
    -o $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_SMD_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB


PATHWAY="UGT"
awk '$3=="RB2011" || $3=="RB2550" || $3=="RB2055" || $3=="UGT49" || $3=="UGT60" || $3=="PD1074" || NR==1' ${DESIGN} > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
    --wide $READY/${VAR}_analysisReady_rank_sbys.tsv \
    --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
    -id uniqueID \
    --study "batch" \
    --treatment "strain" \
    --contrast "RB2011,RB2550,RB2055,UGT49,UGT60,PD1074" \
    --report $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_report.tsv \
    -o $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_SMD_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

PATHWAY="TCA"
awk '$3=="KJ550" || $3=="RB2347" || $3=="VC1265" || $3=="AUM2073" || $3=="VC2524" || $3=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
   --wide $READY/${VAR}_analysisReady_rank_sbys.tsv \
   --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
   -id uniqueID \
   --study "batch" \
   --treatment "strain" \
   --contrast "KJ550,RB2347,VC1265,AUM2073,VC2524,PD1074" \
    --report $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_report.tsv \
    -o $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_SMD_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

PATHWAY="NI"
awk '$3=="CB4856" || $3=="N2" || $3=="CX11314" || $3=="DL238" || $3=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
   --wide $READY/${VAR}_analysisReady_rank_sbys.tsv \
   --design ${ROZ}/dsgn_model_${PATHWAY}.tsv \
   -id uniqueID \
   --study "batch" \
   --treatment "strain" \
   --contrast "CB4856,N2,CX11314,DL238,PD1074" \
    --report $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_report.tsv \
    -o $OUTD/MA_FE_rank_byPath_${VAR}_${PATHWAY}_SMD_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

END
##rm -r $ROZ
