#!/bin/bash
#SBATCH --job-name=sr_kopp
#SBATCH --mail-user=ashley.myrick@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/sr_kopp.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --time=24:00:00
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

module purge
#module load python/3.8

##Create directories
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/sqantiReads_vs_native_kopp_facet_sex
    mkdir -p $OUTD

DESIGN_FILE=$PROJ/design_files/list_samples.csv
DESIGN=$PROJ/design_files/list_samples_ser.csv

echo "sample" > $DESIGN
cat $DESIGN_FILE >> $DESIGN

SQANTI=${PROJ}/sqanti_vs_native/
IDUJC=${PROJ}/

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dser11_2_dser1_ujc.gtf

##Make design files
awk -F, -v sqanti="$SQANTI" -v idujc="$IDUJC" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1;
    print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN" > "$OUTD/kopp_design_4_plots_sex.csv"

##Plots
export PATH=/blue/mcintyre/nkeil/conda_env/sqanti_reads_plots/bin:$PATH

##Run python script


python ${SCRIPTS}/sqanti_reads_testing_13ndk.py \
    -d $OUTD/kopp_design_4_plots_sex.csv \
    -r ${REF_GTF} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p ser_facet_sex1

