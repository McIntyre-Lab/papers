#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cntM_Mo
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/cntM_Mo_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1


### mo17 cau reference
### submit mapping count script (submit_count_mapped) for all samples


## mo17
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/mo17/amb/19_mo17_amb.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/mo17/oz/21-2_mo17_oz.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/mo17/oz/21_mo17_oz.polished.all.hq.mapped.sam

## b73
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/b73/oz/46_b73_oz.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/b73/amb/42_b73_amb.polished.all.hq.mapped.sam

## nc338
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/nc338/amb/89_nc338_amb.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/nc338/oz/96_nc338_oz.polished.all.hq.mapped.sam

## c123
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/c123/amb/113_c123_amb.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/c123/oz/120_c123_oz.polished.all.hq.mapped.sam

## hp301
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/hp301/oz/70_hp301_oz.polished.all.hq.mapped.sam
sh submit_count_mapped_mo17_cau.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_mo17_cau/hp301/amb/67_hp301_amb.polished.all.hq.mapped.sam

cat /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_minimap_mo17_cau_qc_cnt.csv \
    > /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/minimap_mo17_cau_qc_cnt.csv

rm /ufrc/mcintyre/share/maize_ainsworth/pacbio_analysis/counting_recheck/*_minimap_mo17_cau_qc_cnt.csv
