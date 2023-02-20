#!/bin/sh

export PATH=$HOME/conda/envs/mylab/bin:$PATH

## Make flags for chip rna supplemental file

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
IND=$PROJ/RNAseq/Coverage_counts/gene_cvrg_cnts_ortho

## Make output directory
OUTD=$PROJ/figures/mel_sim_ortho_map_scatter
     mkdir -p ${OUTD}

## Prepare data for scatter plots
python $PROJ/scripts/integrated_chip_rna/prep_mel_sim_ortho_map_4_scatter_02avn.py \
    -m ${IND}/cvrGene_ortho2mel_log_uq_apn.csv \
    -s ${IND}/cvrGene_ortho2sim_log_uq_apn.csv \
    -o $PROJ/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv \
    -d ${OUTD}


## Scatter plots of the mel vs sim F/M ratios (avg UQ APN) for
##     conserved male-biased genes and conserved female-biased genes
for COORD in mel sim; do
    python $PROJ/scripts/integrated_chip_rna/plot_mel_sim_ortho_map_scatter.py \
        -m ${OUTD}/ortho_map2${COORD}_conserved_M.csv \
        -f ${OUTD}/ortho_map2${COORD}_conserved_F.csv \
        -p ${OUTD}/ortho_map2${COORD}

    ## Plot the species in COORD but across the mel and sim genome mappings
    python $PROJ/scripts/integrated_chip_rna/plot_mel_sim_ortho_map_scatter.py \
        -m ${OUTD}/ortho_map2bothCoord_conserved_M.csv \
        -s ${COORD} \
        -f ${OUTD}/ortho_map2bothCoord_conserved_F.csv \
        -a ${OUTD}/ortho_map2bothCoord_all.csv \
        -p ${OUTD}/ortho_${COORD}2bothCoord
done
