#!/bin/bash
#SBATCH --job-name=plot_TO
#SBATCH --mail-user=ashley.myrick@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/plot_TO.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --time=06:00:00
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

module purge
#module load python/3.8

##Create directories
PROJ=/blue/mcintyre/share/transcript_ortholog
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/sqantiReads_facet_sex
    mkdir -p $OUTD

DESIGN_FILE=$PROJ/design_files/list_samples.csv

DESIGN_YAK=$PROJ/design_files/list_samples_dyak.csv
DESIGN_SAN=$PROJ/design_files/list_samples_dsan.csv

echo "sample" > "$DESIGN_YAK"
echo "sample" > "$DESIGN_SAN"

awk -v yak="$DESIGN_YAK" -v san="$DESIGN_SAN" '{if ($0 ~/^dyak/) print >> yak; else if ($0 ~/^dsan/) print >> san}' "$DESIGN_FILE"

SQANTI2=${PROJ}/sqanti3QC_yak_san_dorReads_5species_dyak2_dsan1_ujc/
IDUJC2=${PROJ}/

REF_GTF_YAK=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc.gtf
REF_GTF_SAN=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dsan1_ujc.gtf

##Make design files
awk -F, -v sqanti="$SQANTI2" -v idujc="$IDUJC2" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1; print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN_YAK" > "$OUTD/TO_dyak_design_4_plots_comb_TRs_sex.csv"

awk -F, -v sqanti="$SQANTI2" -v idujc="$IDUJC2" 'BEGIN {OFS=","}
NR==1 {if($1=="sample") $1="sampleID"; print $0, "classification_file", "junction_file", "index_file", "sex"}
NR>1 {
    split($1, fields, "_");
    sex=fields[2];
    sample_id=$1; print sample_id, sqanti$1"_classification.txt", sqanti$1"_junctions.txt", idujc$1"_ujc_xscript_link.csv", sex
}' "$DESIGN_SAN" > "$OUTD/TO_dsan_design_4_plots_comb_TRs_sex.csv"


##Plots
export PATH=/blue/mcintyre/nkeil/conda_env/sqanti_reads_plots/bin:$PATH

##Run python script

##yak
python ${SCRIPTS}/sqanti_reads_testing_12ndk.py \
    -d $OUTD/TO_dyak_design_4_plots_comb_TRs_sex.csv \
    -r ${REF_GTF_YAK} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p TO_dyak_comb_TRs_facet_sex

##san
python ${SCRIPTS}/sqanti_reads_testing_12ndk.py \
    -d $OUTD/TO_dsan_design_4_plots_comb_TRs_sex.csv \
    -r ${REF_GTF_SAN} \
    -f sex \
    -o ${OUTD} \
    --all-tables \
    --pca-tables \
    -p TO_dsan_comb_TRs_facet_sex

