#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=flag_rna
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/flag_feature_rna_%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb

##### Make flags for RNA sex specificity #####
PROJ=/ufrc/mcintyre/share/etoh_srna
OUTPUT=$PROJ/rnaseq/sample_level_summary
ROZ=${OUTPUT}/roz_features_rna
SPECIES1="mel"
SPECIES2="sim"
UQFF=$PROJ/rnaseq/uq_ff_map_ff.tsv

mkdir -p ${OUTPUT}

for FEATURE in fragment fusion intron TSS300bpWindow 5UTR 3UTR intergenic; do

    CVR=$PROJ/rnaseq/coverage_cnts/${FEATURE}_coverage_cnts
    echo "#### ${FEATURE} ####
"
    ## Flag on/off
    echo -e "\tFlag on/off..."
    source $PROJ/scripts/rnaseq/run_feature_on_off.sh

    ## Get UQ normalized f/m ratios
    echo -e "\tGet UQ normalized f/m ratios..."
    source $PROJ/scripts/rnaseq/run_feature_f_m_ratio.sh

done

for SPECIES in ${SPECIES1} ${SPECIES2}; do
    ## Combine flags for features of interest
    head -1 ${OUTPUT}/${SPECIES}_fragment_on_off_uq_ff.csv \
        > ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_on_off_uq_ff.csv
    head -1 ${OUTPUT}/${SPECIES}_fusion_on_off_uq_ff.csv \
        > ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_on_off_uq_ff.csv
    head -1 ${OUTPUT}/${SPECIES}_fragment_RNA_sex_bias.csv \
        > ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_RNA_sex_bias.csv
    head -1 ${OUTPUT}/${SPECIES}_fusion_RNA_sex_bias.csv \
        > ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_RNA_sex_bias.csv
    for feature in intron TSS300bpWindow 5UTR 3UTR intergenic; do
        awk 'NR!=1' ${OUTPUT}/${SPECIES}_${feature}_on_off_uq_ff.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_on_off_uq_ff.csv
        awk 'NR!=1' ${OUTPUT}/${SPECIES}_${feature}_on_off_uq_ff.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_on_off_uq_ff.csv
        awk 'NR!=1' ${OUTPUT}/${SPECIES}_${feature}_RNA_sex_bias.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_RNA_sex_bias.csv
        awk 'NR!=1' ${OUTPUT}/${SPECIES}_${feature}_RNA_sex_bias.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_RNA_sex_bias.csv
    done
    awk 'NR!=1' ${OUTPUT}/${SPECIES}_fragment_on_off_uq_ff.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_on_off_uq_ff.csv
    awk 'NR!=1' ${OUTPUT}/${SPECIES}_fusion_on_off_uq_ff.csv \
            >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_on_off_uq_ff.csv
    awk 'NR!=1' ${OUTPUT}/${SPECIES}_fragment_RNA_sex_bias.csv \
        >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_frag_intr_inter_RNA_sex_bias.csv
    awk 'NR!=1' ${OUTPUT}/${SPECIES}_fusion_RNA_sex_bias.csv \
        >> ${OUTPUT}/${SPECIES}_5U_3U_TSS_fus_intr_inter_RNA_sex_bias.csv
done
