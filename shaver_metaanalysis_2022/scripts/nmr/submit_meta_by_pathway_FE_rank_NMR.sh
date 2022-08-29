#!/bash/bin


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPT=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts

DATAIN=$PROJ/nmr
OUTD=$PROJ/nmr/meta_analysis
    mkdir -p $OUTD

for TYPE in cdcl3 d2o
#for TYPE in d2o
do

    DESIGN=$DATAIN/dsgn_${TYPE}_sbys.tsv
        #sampleID        oldSampleID     set     batch   rack_position   run_order       tube_sample     genotype        wormgrowth_sample_name  Solvent parameters      rank    >
        #aos100_ga_ms2_73        aos100_ga_ms2_NMR_CDCL3_73      2       4       C3      73      nmr68   KJ550   aos100_ga_ms2   CDCl3   proton  1       g2      aos100
        #aos101_ga_ms2_121       aos101_ga_ms2_NMR_CDCL3_121     3       6       C3      121     nmr111  UGT49   aos101_ga_ms2   CDCl3   proton  1       g3      aos101

    MODEL=MA_fixedEffect_pathway_NMR_${TYPE}_rank

    WIDE=$DATAIN/ready_${TYPE}_rank_sbys.tsv

    ROZ=$PROJ/roz_MA_nmr_${TYPE}_amm
        mkdir -p $ROZ


## pull out following strain for meta pathway (VC1265=TCA,DL2238=NI,UGT60(VC2512)=UGT
    PATHWAY="ALL"
    awk '$8=="CB4856" || $8=="CX11314" || $8=="DL238" || $8=="VC1265" || $8=="UGT60" || $8=="RB2011" || $8=="RB2550" || $8=="RB2055" || $8=="UGT49" || $8=="KJ550" || $8=="RB2347" || $8=="AUM2073" || $8=="VC2524" || $8=="N2" || $8=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
    python $SCRIPT/meta_analysis.py \
        --wide $WIDE \
        --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
        -id ppm \
        --study "batch" \
        --treatment "genotype" \
        --contrast "AUM2073,CB4856,CX11314,DL238,KJ550,N2,RB2011,RB2055,RB2347,RB2550,UGT49,UGT60,VC1265,VC2524,PD1074"  \
        --forest no \
        --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
        -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
        --model FE \
        --effectSize SMD \
        --cmMethod UB

: <<'END'

## pull out following strain for meta pathway (VC1265=TCA,DL2238=NI,UGT60(VC2512)=UGT
    PATHWAY="NInoN2"
    awk '$8=="CB4856" || $8=="CX11314" || $8=="DL238" || $8=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
    python $SCRIPT/meta_analysis.py \
        --wide $WIDE \
        --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
        -id ppm \
        --study "batch" \
        --treatment "genotype" \
        --contrast "CB4856,CX11314,DL238,PD1074" \
        --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
        -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
        --model FE \
        --effectSize SMD \
        --cmMethod UB

## pull out following strain for meta pathway (VC1265=TCA,DL2238=NI,UGT60(VC2512)=UGT
PATHWAY="COMBO"
awk '$8=="VC1265" || $8=="DL238" || $8=="UGT60" || $8=="PD1074" || NR==1' ${DESIGN} > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
python $SCRIPT/meta_analysis.py \
    --wide $WIDE \
    --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
    -id ppm \
    --study "batch" \
    --treatment "genotype" \
    --contrast "VC1265,DL238,UGT60,PD1074" \
    --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
    -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
    --model FE \
    --effectSize SMD \
    --cmMethod UB

    ## pull out strains with mutationsin UGT pathway
    PATHWAY="UGT"
    awk '$8=="RB2011" || $8=="RB2550" || $8=="RB2055" || $8=="UGT49" || $8=="UGT60" || $8=="PD1074" || NR==1' ${DESIGN} > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
    python $SCRIPT/meta_analysis.py \
        --wide $WIDE \
        --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
        -id ppm \
        --study "batch" \
        --treatment "genotype" \
        --contrast "RB2011,RB2550,RB2055,UGT49,UGT60,PD1074" \
        --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
        -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
        --model FE \
        --effectSize SMD \
        --cmMethod UB

    PATHWAY="TCA"
    awk '$8=="KJ550" || $8=="RB2347" || $8=="VC1265" || $8=="AUM2073" || $8=="VC2524" || $8=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
    python $SCRIPT/meta_analysis.py \
        --wide $WIDE \
        --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
        -id ppm \
        --study "batch" \
        --treatment "genotype" \
        --contrast "KJ550,RB2347,VC1265,AUM2073,VC2524,PD1074" \
        --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
        -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
        --model FE \
        --effectSize SMD \
        --cmMethod UB

    PATHWAY="NI"
    awk '$8=="CB4856" || $8=="N2" || $8=="CX11314" || $8=="DL238" || $8=="PD1074" || NR==1' $DESIGN > ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv
    python $SCRIPT/meta_analysis.py \
        --wide $WIDE \
        --design ${ROZ}/dsgn_model_${TYPE}_${PATHWAY}.tsv \
        -id ppm \
        --study "batch" \
        --treatment "genotype" \
        --contrast "CB4856,N2,CX11314,DL238,PD1074" \
        --report $OUTD/${MODEL}_${PATHWAY}_report.tsv \
        -o $OUTD/${MODEL}_${PATHWAY}_summary.tsv \
        --model FE \
        --effectSize SMD \
        --cmMethod UB
END
done

#rm -r $ROZ
