#!/bash/bin

## lactic acid forest plot for d2o data (note lactic acid NOT in cdcl3 dataset!!!

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPT=/TB14/TB14/github/SECIMTools/src/scripts

DATAIN=$PROJ/nmr

OUTD=$PROJ/nmr/meta_analysis
    mkdir -p $OUTD

DESIGN=$DATAIN/dsgn_d2o_sbys.tsv
        #sampleID        oldSampleID     set     batch   rack_position   run_order       tube_sample     genotype        wormgrowth_sample_name  Solvent parameters      rank    >
        #aos100_ga_ms2_73        aos100_ga_ms2_NMR_CDCL3_73      2       4       C3      73      nmr68   KJ550   aos100_ga_ms2   CDCl3   proton  1       g2      aos100
        #aos101_ga_ms2_121       aos101_ga_ms2_NMR_CDCL3_121     3       6       C3      121     nmr111  UGT49   aos101_ga_ms2   CDCl3   proton  1       g3      aos101

WIDE=$DATAIN/ready_d2o_rank_sbys.tsv

ROZ=$PROJ/roz_MA_nmr_d2o_amm
    mkdir -p $ROZ

export PATH=/TB14/TB14/galaxy_18feb21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH

## subset data to only ppms sig in meta pathway and not in individual NI strains (omitted N2)
head -n1 $DATAIN/ready_d2o_rank_sbys.tsv > $ROZ/plots_d2o_rank_sbys.tsv
grep -E '(ppm_1_724|ppm_1_8979|ppm_1_9349|ppm_2_0967|ppm_2_215|ppm_2_3291|ppm_2_7514|ppm_3_029|ppm_3_2672|ppm_3_8184|ppm_3_9741|ppm_4_252|ppm_4_357|ppm_4_5912|ppm_4_6416|ppm_4_6629|ppm_5_5153|ppm_5_5317|ppm_5_7596|ppm_5_8514|ppm_5_8634|ppm_5_8836|ppm_6_703|ppm_6_8322|ppm_7_2579|ppm_7_5467|ppm_7_6519|ppm_7_7129|ppm_8_042)' $DATAIN/ready_d2o_rank_sbys.tsv >> $ROZ/plots_d2o_rank_sbys.tsv

## meta pathway 
PATHWAY="NI_no_N2"
awk '$8=="CB4856" || $8=="CX11314" || $8=="DL238" ||  $8=="PD1074" || NR==1' ${DESIGN} > ${ROZ}/dsgn_model_d2o_NI_no_N2.tsv
python $SCRIPT/meta_analysis.py \
    --wide $ROZ/plots_d2o_rank_sbys.tsv \
    --design ${ROZ}/dsgn_model_d2o_NI_no_N2.tsv \
    -id ppm \
    --study "batch" \
    --treatment "genotype" \
    --contrast "CB4856,CX11314,DL238,PD1074" \
    --report $OUTD/MA_FE_path_d20_NInoN2_sigMP_notSigMS_report.tsv \
    -o $OUTD/MA_FE_path_d2o_NInoN2_sigMP_notSigMS_summary.tsv \
    --forest $OUTD/MA_FE_path_d20_NInoN2_sigMP_notSigMS \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

rm -r $ROZ
