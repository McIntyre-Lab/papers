#!/bin/bash
#SBATCH --job-name=catujcTO
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/catujcTO_%A-%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2gb
#SBATCH --time=4:00:00
#SBATCH --account=mcintyre
##SBATCH --array=1
#SBATCH --array=1-2

module purge

ml conda
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

PROJ=/blue/mcintyre/share/transcript_ortholog

IND=$PROJ/ujc_byGene_output/
OUTD=$PROJ

## Get design file information (2 rows, dsan and dyak)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

#ROZ=$PROJ/roz_${SPECIES}
#    mkdir -p $ROZ

echo "slurm temp is here:  $SLURM_TMPDIR"

echo "Catting xscript link"
# cat xscript link
echo "source,transcriptID,geneID,jxnHash,jxnString" > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_simpleCat.csv
tail -n +2 -q $IND/${SPECIES}_data_*_ujc_xscript_link.csv >> $OUTD/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_simpleCat.csv

# echo "Catting count files"
# cat count files
# echo "source,jxnHash,numTranscripts" > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_count_simpleCat.csv
# tail -n +2 -q $IND/${SPECIES}_data_*_ujc_count.csv >> $OUTD/${SPECIES}_data_2_${GENOME}_ujc_count_simpleCat.csv

echo "Catting dsc files"
# cat descript
echo "jxnHash,geneID,flagMultiXscript,flagMultiGene,numJxn,chr,strand,donorStart,acceptorEnd,jxnString" > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_dscrptn_simpleCat.csv
tail -n +2 -q $IND/${SPECIES}_data_*_ujc_dscrptn.csv >> $OUTD/${SPECIES}_data_2_${GENOME}_ujc_dscrptn_simpleCat.csv

echo "Catting gtf"
## cat gtf (no header!)
cat $IND/${SPECIES}_data_*_ujc.gtf > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_simpleCat.gtf

#echo "Sorting files..."
## SORT xscript
#(head -n 1 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_xscript_link.csv && tail -n +2 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_xscript_link.csv | sort -t ',' -k1,1 -k3,3) \
#    > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_xscript_link.csv
#
## sort by gene then jxnHash
#(head -n 1 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv && tail -n +2 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv | sort -t ',' -k2,2 -k1,1) \
#    > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv
#
## sort by sample
#(head -n 1 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_count.csv && tail -n +2 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc_count.csv | sort -t ',' -k1,1) \
#    > $OUTD/${SPECIES}_data_2_${GENOME}_ujc_count.csv
#
## sort gtf by chrom, start, jnxHash (transcriptID)
#sort -t $'\t' -k1,1 -k4,4n -k10,10 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_ujc.gtf > $OUTD/${SPECIES}_data_2_${GENOME}_ujc.gtf

echo "Completed~!"
