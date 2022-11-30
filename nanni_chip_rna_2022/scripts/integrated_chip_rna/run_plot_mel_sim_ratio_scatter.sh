#!/bin/sh

export PATH=$HOME/conda/envs/mylab/bin:$PATH

## Make flags for chip rna supplemental file

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
IND=$PROJ/supp_files

## Make output directory
OUTD=$PROJ/figures/mel_sim_ratio_scatter
     mkdir -p ${OUTD}

## Prepare data for scatter plots of conserved sex bias
python $PROJ/scripts/integrated_chip_rna/prep_mel_sim_ratio_4_scatter.py \
    -mf $PROJ/RNAseq/model_output/mel_frag_flags_kitchen_sink.csv \
    -ma ~/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_exon_fragment_annotations.csv \
    -sf $PROJ/RNAseq/model_output/sim_frag_flags_kitchen_sink.csv \
    -sa ~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/flybase202/event_analysis_annotations/dsim202_annotations_150bp_reads/dsim202_150bp_exon_fragment_annotations.csv \
    -o $PROJ/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv \
    -d ${OUTD}

## Scatter plots of the mel vs sim F/M ratios (avg UQ APN) for
##     conserved male-biased genes and conserved female-biased genes
python $PROJ/scripts/integrated_chip_rna/plot_mel_sim_ratio_scatter_04avn.py \
    -m ${OUTD}/ortho_MF_ratio_conserved_M.csv \
    -f ${OUTD}/ortho_MF_ratio_conserved_F.csv \
    --Mmel-Fsim ${OUTD}/ortho_MF_ratio_M_mel_F_sim.csv \
    --Fmel-Msim ${OUTD}/ortho_MF_ratio_F_mel_M_sim.csv \
    -d ${OUTD}

## Prepare data for scatter plots of mel and sim ratios of all genes
    # Top right quadrant will be conserved male (1 - F/M, both species)
    # Bottom left quadrant will be conserved female (M/F - 1, both species)
    # Bottom right quadrant will be mel male and simulans female
    # Bottom left quadrant will be sim male and mel female
python $PROJ/scripts/integrated_chip_rna/prep_mel_sim_ratio_4_scatter_allGenes.py \
    -mf $PROJ/RNAseq/model_output/mel_frag_flags_kitchen_sink.csv \
    -ma ~/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_exon_fragment_annotations.csv \
    -sf $PROJ/RNAseq/model_output/sim_frag_flags_kitchen_sink.csv \
    -sa ~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/flybase202/event_analysis_annotations/dsim202_annotations_150bp_reads/dsim202_150bp_exon_fragment_annotations.csv \
    -o $PROJ/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv \
    -d ${OUTD}

## Plot scatter plot of all gene ratios (unbiased included)
python $PROJ/scripts/integrated_chip_rna/plot_mel_sim_ratio_scatter_allGenes_02avn.py \
    -i ${OUTD}/ortho_MF_ratio_all_gene.csv \
    -d ${OUTD}
