
#!/bin/bash

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH


PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

READY=$PROJ/SLAW_UGA_Output/analysisReady_datasets


DESIGN=$PROJ/design_files/dsgn_GT_${VARDIR}_pairs_slaw.tsv
    #sampleID        sample  strain  batch   sampleType      newBatch        sampleType_batch        strain_w_poolPD pair    set
    #b1_aos_122      aos_122 RB2011  1       strain  b1      strain_b1       RB2011  1       set1
    #b1_aos_127      aos_127 CX11314 1       strain  b1      strain_b1       CX11314 1       set1

OUTD=$PROJ/SLAW_UGA_Output/meta_analysis_${VAR}_rerun
    mkdir -p $OUTD

ROZ=$PROJ/SLAW_UGA_Output/roz_MA/${VAR}
    mkdir -p $ROZ


for set in set1 set2 set3
do
    # for each mutant in setX subset 
    for mut in $(awk -v set=${set} '$2==set {print $1}' $PROJ/design_files/batch_model_tests/all_mutantsList.tsv)
    do
        ## pull out mutant and all pd1074 in the set from DESIGN, output to new design file for input into model
        awk -F "\t"  -v mut=$mut -v set=$set '($3==mut || $3=="PD1074" || NR==1) && ($10==set || NR==1)' ${DESIGN} > ${ROZ}/dsgn_${set}_${mut}.tsv

	python $SCRIPTS/meta_analysis.py \
	    --wide $READY/${VAR}_analysisReady_rank_sbys.tsv \
	    --design ${ROZ}/dsgn_${set}_${mut}.tsv \
	    -id uniqueID \
	    --study "batch" \
	    --treatment "strain" \
	    --contrast "$mut,PD1074" \
            --forest no \
            --report $OUTD/MA_FE_rank_byMut_${VAR}_${mut}_report.txt \
	    -o $OUTD/MA_FE_rank_byMut_${VAR}_${mut}_SMD_summary.tsv \
            --model FE \
            --effectSize SMD \
            --cmMethod UB

    done
done

