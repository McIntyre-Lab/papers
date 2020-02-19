#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=df_1Line
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/df_1Line%A_%a.out
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1


## split design file by sampleID

PROJ=/ufrc/mcintyre/share/maize_ainsworth

DESIGN=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
OUT=$PROJ/design_files


for ID in 19 113 46 89 21 70 96 120 67 42 21-2
do
    grep ",$ID," $DESIGN > $OUT/df_maize_test_PacBio_path_noHeader_${ID}.csv
    ## add in full path
    sed 's/PacBio/\/ufrc\/mcintyre\/share\/maize_ainsworth\/original_data\/PacBio/g' $OUT/df_maize_test_PacBio_path_noHeader_${ID}.csv > $OUT/df_maize_test_PacBio_fullpath_noHeader_${ID}.csv
done
