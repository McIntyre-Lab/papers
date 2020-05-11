#!/bin/bash


## for isolates snps called against CC group ref:
	## remove header for import into sas


### Set Directories
PROJ=/home/ammorse/TB14/staph_relapse
VCF_IN=${PROJ}/vcf_w_targets_file

#create list to loop over
cd ${VCF_IN}
ls *_snps.vcf | sed 's/.vcf//g' > $PROJ/design_files/list_vcf_4_sas.txt

while IFS= read -r VCF
do

    sed '/^##/d' ${VCF_IN}/${VCF}.vcf > ${VCF_IN}/${VCF}_noHeader.vcf

done < $PROJ/design_files/list_vcf_4_sas.txt
