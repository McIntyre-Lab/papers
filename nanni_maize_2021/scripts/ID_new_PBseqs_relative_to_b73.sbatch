#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=extract
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/extract_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6M
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --array=1

# for samples NOT b73, identify transcripts that do NOT map to ref B73 but DO map to ref Mo17
#	these are following rows in DF: 1,3,-9,11


### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
INPUT=$PROJ/uncollapsed_b73_IsoSeq_reads

OUTPUT=$PROJ/evidence_genotype_specificity
    mkdir -p $OUTPUT

## list mo17 transcripts that do NOT map to b73
#	list_19_mo17_amb.ignored_unmapped_ids.b73_ref.txt
#	list_21-2_mo17_oz.ignored_unmapped_ids.b73n_ref.txt
#	list_21_mo17_oz.ignored_unmapped_ids.b73_ref.txt

## do transcripts in the above lists map?
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/amb/19_mo17_amb.polished.all.hq.fasta
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/oz/21-2_mo17_oz.polished.all.hq.fasta
#	should exist in $PROJ/mapping_minimap2_mo17_yan/mo17/oz/21_mo17_oz.polished.all.hq.fasta

## Design file
DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE=${ID}_${GENO}_${TRT}

LIST=$INPUT/list_${SAMPLE}.ignored_unmapped_ids.b73_ref.txt
FASTA=$PROJ/mapping_minimap2_mo17_yan/${GENO}/${TRT}/${SAMPLE}.polished.all.hq.fasta

echo -e " list of unmapped to B73 to check and see if maps to Mo17:  $LIST \n
          fasta file to search:  $FASTA \n"

while read ID
do
    echo -e "looking for ID:  $ID\n"
    grep $ID $FASTA >> $OUTPUT/${SAMPLE}_new_relative_to_B73.txt 
done 
