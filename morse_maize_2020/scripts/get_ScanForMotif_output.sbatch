#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=scanformotif
#SBATCH --output=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB/scripts/SLURM_LOGS/scanformotif_%A.out
#SBATCH --qos=mcintyre
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1

module purge

## Get ScanForMotif webserver output

### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
OUTPUT=$PROJ/scanformotif
    mkdir -p ${OUTPUT}

## Set species
SPECIES=maize

# 00
#wget -O ${OUTPUT}/${SPECIES}.00.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.44.08_2019_47782_report.txt

# 01
#wget -O ${OUTPUT}/${SPECIES}.01.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.44.23_2019_73989_report.txt

# 02
#wget -O ${OUTPUT}/${SPECIES}.02.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.44.31_2019_30759_report.txt

# 03
wget -O ${OUTPUT}/${SPECIES}.03.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.44.57_2019_55336_report.txt

# 04
wget -O ${OUTPUT}/${SPECIES}.04.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.45.11_2019_66033_report.txt

# 05
wget -O ${OUTPUT}/${SPECIES}.05.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.45.23_2019_63245_report.txt

# 06
wget -O ${OUTPUT}/${SPECIES}.06.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.48.46_2019_03980_report.txt

# 07
wget -O ${OUTPUT}/${SPECIES}.07.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.48.53_2019_28303_report.txt

# 08
wget -O ${OUTPUT}/${SPECIES}.08.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.49.04_2019_61652_report.txt

# 09
wget -O ${OUTPUT}/${SPECIES}.09.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.49.08_2019_37822_report.txt

# 10
wget -O ${OUTPUT}/${SPECIES}.10.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.50.27_2019_44916_report.txt

# 11
wget -O ${OUTPUT}/${SPECIES}.11.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.51.02_2019_37833_report.txt

# 12
wget -O ${OUTPUT}/${SPECIES}.12.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.52.45_2019_59820_report.txt

# 13
wget -O ${OUTPUT}/${SPECIES}.13.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.54.00_2019_53000_report.txt

# 14
wget -O ${OUTPUT}/${SPECIES}.14.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.54.22_2019_46001_report.txt

# 15
wget -O ${OUTPUT}/${SPECIES}.15.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.54.22_2019_26030_report.txt

# 16
wget -O ${OUTPUT}/${SPECIES}.16.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.55.28_2019_55175_report.txt

# 17
wget -O ${OUTPUT}/${SPECIES}.17.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.55.47_2019_53609_report.txt

# 18
wget -O ${OUTPUT}/${SPECIES}.18.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_16.14.21_2019_20000_report.txt

# 19
wget -O ${OUTPUT}/${SPECIES}.19.SFM.out http://bioanalysis.otago.ac.nz/sfm/tmp/Fri_Oct__4_15.58.36_2019_77081_report.txt
