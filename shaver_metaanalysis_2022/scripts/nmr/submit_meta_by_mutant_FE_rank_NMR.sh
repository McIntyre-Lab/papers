
#!/bin/bash

###MODEL=MA_fixedEffect_NMR_cdcl3_rank

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPT=/TB14/TB14/github/SECIMTools/src/scripts

DATAIN=$PROJ/nmr

OUTD=$PROJ/nmr/meta_analysis
    mkdir -p $OUTD

for TYPE in cdcl3 d2o
do

    DESIGN=$DATAIN/dsgn_${TYPE}_sbys.tsv
	#sampleID        oldSampleID     set     batch   rack_position   run_order       tube_sample     genotype        wormgrowth_sample_name  Solvent parameters      rank    study_group     sample
	#aos100_ga_ms2_73        aos100_ga_ms2_NMR_CDCL3_73      2       4       C3      73      nmr68   KJ550   aos100_ga_ms2   CDCl3   proton  1       g2      aos100
	#aos101_ga_ms2_121       aos101_ga_ms2_NMR_CDCL3_121     3       6       C3      121     nmr111  UGT49   aos101_ga_ms2   CDCl3   proton  1       g3      aos101

    MODEL=MA_fixedEffect_NMR_${TYPE}_rank

    WIDE=$DATAIN/ready_${TYPE}_rank_sbys.tsv

    ROZ=$PROJ/roz_MA_nmr_${TYPE}_amm
        mkdir -p $ROZ

    export PATH=/TB14/TB14/galaxy_18feb21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH

    for set in 1 2 3
    do
        # for each mutant in setX subset 
        for mut in $(awk -v set=set${set} '$2==set {print $1}' $PROJ/design_files/batch_model_tests/all_mutantsList.tsv) 
        do
            echo "mutant is $mut and set is ${set}"
            ## pull out mutant and all pd1074 in the set from DESIGN, output to new design file for input into model
            awk -v mut=$mut -v set="$set" '($8==mut || $8=="PD1074" || NR==1) && ($3==set || NR==1)' ${DESIGN} > ${ROZ}/dsgn_nmr_${TYPE}_set${set}_${mut}.tsv

    	    python $SCRIPT/meta_analysis.py \
	        --wide $WIDE \
    	        --design ${ROZ}/dsgn_nmr_${TYPE}_set${set}_${mut}.tsv \
	        -id ppm \
                --study "batch" \
	        --treatment "genotype" \
	        --contrast "$mut,PD1074" \
                --report $OUTD/${MODEL}_${mut}_${TYPE}_report.txt \
	        -o $OUTD/${MODEL}_${mut}_${TYPE}_SMD_summary.tsv \
                --model FE \
                --effectSize SMD \
                --cmMethod UB
        done
    done
done
