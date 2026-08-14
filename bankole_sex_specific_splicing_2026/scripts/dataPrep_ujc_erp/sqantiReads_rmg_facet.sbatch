#!/bin/bash
#SBATCH --job-name=sr_rmg
#SBATCH --mail-user=ashley.myrick@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/sr_rmg.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --time=24:00:00
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

module purge
#module load python/3.8

##Create directories
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/sqantiReads_facet_sex
    mkdir -p $OUTD

DESIGN_FILE=$PROJ/design_files/list_samples.csv

DESIGN_MEL=$PROJ/design_files/list_samples_mel.csv
DESIGN_SIM=$PROJ/design_files/list_samples_sim.csv


SQANTI=${PROJ}/sqanti3QC_ser_dorReads_2_5species/
IDUJC=${PROJ}/

REF_GTF_MEL=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf
REF_GTF_SIM=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dsim2_ujc.gtf

##Make design files
awk -F, -v sqanti="$SQANTI" -v idujc="$IDUJC" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1;
    print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN_MEL" > "$OUTD/rmg_mel_design_4_plots_comb_TRs_sex.csv"

awk -F, -v sqanti="$SQANTI" -v idujc="$IDUJC" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1;
    print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN_SIM" > "$OUTD/rmg_sim_design_4_plots_comb_TRs_sex.csv"


##Plots
export PATH=/blue/mcintyre/nkeil/conda_env/sqanti_reads_plots/bin:$PATH

##Run python script

##mel
python ${SCRIPTS}/sqanti_reads_testing_13ndk.py \
    -d $OUTD/rmg_mel_design_4_plots_comb_TRs_sex.csv \
    -r ${REF_GTF_MEL} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p rmg_mel_comb_TRs_facet_sex

##sim
python ${SCRIPTS}/sqanti_reads_testing_13ndk.py \
    -d $OUTD/rmg_sim_design_4_plots_comb_TRs_sex.csv \
    -r ${REF_GTF_SIM} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p rmg_sim_comb_TRs_facet_sex

