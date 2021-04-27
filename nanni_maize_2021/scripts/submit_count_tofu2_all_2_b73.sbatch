#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cnt_clus
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/cnt_clus_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1


### submit tofu2 count script (submit_count_tofu2.sh) for all sampleIDs against b73 ref

## mo17
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/19_mo17_amb/19_mo17_amb
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/21_mo17_oz/21_mo17_oz
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/21-2_mo17_oz/21-2_mo17_oz


## b73
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/46_b73_oz/46_b73_oz
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/42_b73_amb/42_b73_amb

## nc338
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/89_nc338_amb/89_nc338_amb
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/96_nc338_oz/96_nc338_oz

## c123
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/113_c123_amb/113_c123_amb
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/120_c123_oz/120_c123_oz

## hp301
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/70_hp301_oz/70_hp301_oz
sh submit_count_tofu2_b73.sh /ufrc/mcintyre/share/maize_ainsworth/tofu2_b73/67_hp301_amb/67_hp301_amb

cat /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_tofu2_b73_counts_cluster_report.txt \
    > /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/tofu2_b73_counts_cluster_report.txt

cat /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_tofu2_b73_counts_collaped_group.txt \
    > /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/tofu2_b73_counts_collaped_group.txt

rm /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_tofu2_b73_counts_cluster_report.txt
rm /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_tofu2_b73_counts_collaped_group.txt
