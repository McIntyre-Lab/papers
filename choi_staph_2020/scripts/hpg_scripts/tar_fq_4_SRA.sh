#!/bin/bash 


## tar staph relapse files together for sra submission
	## note NOT including failed library

PROJ=/ufrc/mcintyre/share/staph_relapse
DF=$PROJ/design_files

cd $PROJ/original_data

tar -v -c -z -f $PROJ/morse_staph_relapse_sra.tgz --files-from=$DF/design_sra_pack_fq.csv 
