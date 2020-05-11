#!/bin/bash


### merge vcf files by reference
### create a list of merged vcf files
### call snps for each isolate using this list of snps

### Set Directories
PROJ=/home/ammorse/TB14/staph_relapse
SCRIPTS=${PROJ}/scripts

OUTPUT=${PROJ}/vcf_by_ref
    mkdir -p ${OUTPUT}


#while IFS="," read -r PAIR1 PAIR2
#do

    ## zip and index starting vcf files
#    bgzip ${PROJ}/freebayes/${PAIR1}_2_CCRef_*_renamed_snps.recode.vcf
#    bgzip ${PROJ}/freebayes/${PAIR2}_2_CCRef_*_renamed_snps.recode.vcf

#    tabix -f -p vcf ${PROJ}/freebayes/${PAIR1}_2_CCRef_*_renamed_snps.recode.vcf.gz
#    tabix -f -p vcf ${PROJ}/freebayes/${PAIR2}_2_CCRef_*_renamed_snps.recode.vcf.gz

#done < ${PROJ}/design_files/isolate_pairs_only_noS15_noHeader.csv


## merge pairs by reference
for REF in ED98 Newman TCH60 ST20130941 CA_347 MSSA476
do
    bcftools merge -m both -o ${OUTPUT}/${REF}_snps.vcf.gz -O z ${PROJ}/freebayes/*_${REF}_renamed_snps.recode.vcf.gz

    bgzip -d ${OUTPUT}/${REF}_snps.vcf.gz

    ## remove header
    sed '/^##/d' ${OUTPUT}/${REF}_snps.vcf > ${OUTPUT}/${REF}_snps_noHeader.vcf

done

## create list of merged vcf files by ref
ls ${OUTPUT}/*_snps_noHeader.vcf > $PROJ/design_files/list_vcf_by_ref.txt

