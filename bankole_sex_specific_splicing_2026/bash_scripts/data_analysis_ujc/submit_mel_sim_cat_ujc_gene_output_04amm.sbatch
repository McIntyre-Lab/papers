#!/bin/bash
#SBATCH --job-name=catujcR
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/catujcR_%A-%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6gb
#SBATCH --time=8:00:00
#SBATCH --account=mcintyre
##SBATCH --array=1
#SBATCH --array=1-2

module purge

ml conda
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

PROJ=/blue/mcintyre/share/rmg_lmm_dros_data

IND=$PROJ/ujc_byGene_output
OUTD=$PROJ

## Get design file information (2 rows, dmel and dsim)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=d${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

echo "Catting xscript link"
# cat xscript link
echo "source,transcriptID,geneID,jxnHash,jxnString" > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_xscript_link.csv
tail -n +2 -q $IND/${SPECIES}_data_*_ujc_xscript_link.csv >> $OUTD/${SPECIES}_data_2_${GENOME}_ujc_xscript_link.csv

echo "Catting dsc files"
# cat descript
echo "jxnHash,geneID,flagMultiXscript,flagMultiGene,numJxn,chr,strand,donorStart,acceptorEnd,jxnString" > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv
tail -n +2 -q $IND/${SPECIES}_data_*_ujc_dscrptn.csv >> $OUTD/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv

echo "Catting gtf"
## cat gtf (no header!)
cat $IND/${SPECIES}_data_*_ujc.gtf > $OUTD/${SPECIES}_data_2_${GENOME}_ujc.gtf

echo "Completed~!"
