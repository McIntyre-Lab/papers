#!/bin/bash


PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing
IN=$PROJ/submission/supplementary_files/fiveSpecies_annotations
OUT=$PROJ/anno_testing
    mkdir -p $OUT
    
export PATH=/nfshome/ammorse/conda/envs/mcintyre_updated/bin:$PATH

x=dmel650
y=dmel6

python $PROJ/scripts/create_fiveSpecies_full_anno_dmel6_file_02amm.py \
    --jxnhash_fbtr_file $IN/fiveSpecies_2_${y}_anno_files/${x}_2_${y}_ujc_xscript_link.csv \
    --genome ${y} \
    --erp_file $IN/fiveSpecies_2_${y}_anno_files/fiveSpecies_2_${y}_ujc_er_vs_fiveSpecies_2_${y}_ujc_infoERP.csv \
    --esp_file $IN/fiveSpecies_2_${y}_anno_files/fiveSpecies_2_${y}_ujc_es_vs_fiveSpecies_2_${y}_ujc_infoESP.csv \
    --flag_file $IN/fiveSpecies_2_${y}_anno_files/flag_fiveSpecies_2_${y}_ujc_w_IR_flag.csv \
    --output_file $OUT/fiveSpecies_${y}_full_annotation.csv \
    --mismatch $OUT/fiveSpecies_${y}_full_geneID_mismatch.csv

x=dsim202
x=dser11
x=dyak21
x=dsan11

y=dsim2
y=dser1
y=dyak2
y=dsan1

# Define the species pairs
x_species=("dsim202" "dser11" "dyak21" "dsan11")
y_species=("dsim2" "dser1" "dyak2" "dsan1")

# Loop over pairs of species
for i in "${!x_species[@]}"
do
    x="${x_species[$i]}"
    y="${y_species[$i]}"

echo "running ${y}"

python $PROJ/scripts/create_fiveSpecies_full_anno_nonMel_files_04amm.py \
    --jxnhash_fbtr_file $IN/fiveSpecies_2_${y}_anno_files/${x}_2_${y}_ujc_xscript_link.csv \
    --genome ${y} \
    --erp_file $IN/fiveSpecies_2_${y}_anno_files/fiveSpecies_2_${y}_ujc_er_vs_fiveSpecies_2_${y}_ujc_infoERP.csv \
    --esp_file $IN/fiveSpecies_2_${y}_anno_files/fiveSpecies_2_${y}_ujc_es_vs_fiveSpecies_2_${y}_ujc_infoESP.csv \
    --flag_file $IN/fiveSpecies_2_${y}_anno_files/flag_fiveSpecies_2_${y}_ujc_w_IR_flag.csv \
    --cross_file $PROJ/cross_species_link_files/fivespecies_${y}_w_geneid_02amm.csv \
    --output_file $OUT/fiveSpecies_${y}_full_annotation.csv \
    --mismatch $OUT/fiveSpecies_${y}_full_geneID_mismatch.csv
done

