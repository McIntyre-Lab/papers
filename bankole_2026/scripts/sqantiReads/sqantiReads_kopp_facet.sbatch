#!/bin/bash
#SBATCH --job-name=plot_kopp
#SBATCH --mail-user=ashley.myrick@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/plot_kopp.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --time=06:00:00
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

module purge
#module load python/3.8

##Create directories
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/sqantiReads_facet_sex
    mkdir -p $OUTD

DESIGN_FILE=$PROJ/design_files/list_samples.csv
DESIGN=$PROJ/design_files/list_samples_ser.csv

echo "sample" > $DESIGN
cat $DESIGN_FILE >> $DESIGN

SQANTI2=${PROJ}/sqanti3QC_ser_dorReads_2_5species_dser1_ujc/
IDUJC2=${PROJ}/

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dser1_ujc.gtf

##Make design files
awk -F, -v sqanti="$SQANTI2" -v idujc="$IDUJC2" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1;
    print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN" > "$OUTD/kopp_dser_design_4_plots_comb_TRs_sex.csv"


##Plots
export PATH=/blue/mcintyre/nkeil/conda_env/sqanti_reads_plots/bin:$PATH

##Run python script

##ser
python ${SCRIPTS}/sqanti_reads_testing_12ndk.py \
    -d $OUTD/kopp_dser_design_4_plots_comb_TRs_sex.csv \
    -r ${REF_GTF} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p kopp_dser_comb_TRs_facet_sex


##Remove temporary design file
rm $DESIGN
