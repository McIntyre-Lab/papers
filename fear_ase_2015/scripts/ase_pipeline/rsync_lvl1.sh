#!/bin/bash
CURRDIR=$HOME/storage/cegs_ase_paper/ase_lvl1_filtered_vcf_files
for FILE in $CURRDIR/*.vcf;
do
    if [ -e $FILE ]
    then
        $HOME/opt/bin/bgzip $FILE
        $HOME/opt/bin/tabix -p vcf ${FILE}.gz
    fi
done

rsync -av $CURRDIR hic:/scratch/lfs/mcintyre/cegs_ase_paper/ase_pipeline_output/
