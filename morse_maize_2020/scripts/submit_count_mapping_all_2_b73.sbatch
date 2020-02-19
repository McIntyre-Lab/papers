#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cnt_map
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/cnt_map_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1


### b73 reference
### submit mapping count script (submit_count_mapped) for all samples

#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/mo17/amb/19_mo17_amb.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/b73/oz/46_b73_oz.polished.all.hq.mapped.sam
sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/nc338/amb/89_nc338_amb.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/c123/amb/113_c123_amb.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/mo17/oz/21_mo17_oz.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/hp301/oz/70_hp301_oz.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/nc338/oz/96_nc338_oz.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/c123/oz/120_c123_oz.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/hp301/amb/67_hp301_amb.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/b73/amb/42_b73_amb.polished.all.hq.mapped.sam
#sh submit_count_mapped_b73.sh /ufrc/mcintyre/share/maize_ainsworth/mapping_minimap2_b73/mo17/oz/21-2_mo17_oz.polished.all.hq.mapped.sam

