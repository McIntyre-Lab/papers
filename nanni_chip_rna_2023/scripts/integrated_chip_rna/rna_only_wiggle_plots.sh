#!/bin/sh

## Get RNA wiggle plots only for genes of interest
## Use text data and R scripts from Jeremy's gene plots

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
INPUT=${PROJ}/giant_manuscript/figures/gene_plots/text_data
OUTPUT=${PROJ}/CHIP_RNA_ms/figures/rna_wiggle_plots
ROZ=${OUTPUT}/roz_rna_wig_files
mkdir -p ${ROZ}

## Genes of interest for mel (sim ortholog)
## dsx FBgn0000504 (GD19448 FBgn0190947)
## fln FBgn0005633 (GD12234 FBgn0183970)
## Act88F FBgn0000047 (Act88F FBgn0012823)
## fru FBgn0004652 (fru FBgn0043395)
## msl-1 FBgn0005617 (GD27005 FBgn0268295)
## msl-2 FBgn0005616 (GD15895 FBgn0187546)
## msl-3 FBgn0002775 (GD13090 FBgn0184812)
## mof FBgn0014340 (GD28183 FBgn0269473)
## mle FBgn0002774 (GD10344 FBgn0182116)
## roX1 FBgn0019661 X:3857230-3862697(-) (unannotated Scf_X:3450629-3454544(-))
## roX2 FBgn0019660 X:11579065-11580432(+) (GD29387 FBgn0270677 Scf_X:10878530-10879855(+))
## ***roX gene names were set as regions so Jeremy's text data identifies them as such:
##     mel roX1 = X_3811070_3872580
##     sim roX1 = Scf_X_3403220_3464520
##     mel roX2 = X_11576780_11585890
##     sim roX2 = Scf_X_10876200_10885400

## Manually determining window based on looking at previous version of the plots

## Change directory for R script
cd ${PROJ}/scripts/integrated_chip_rna
#cd ${PROJ}/scripts

## mel wiggle plots first
SPECIES=mel
for GENE in FBgn0000504 FBgn0005633 FBgn0000047 FBgn0004652 FBgn0005617 FBgn0005616 FBgn0002774 FBgn0014340 FBgn0002775 FBgn0019661 X_11576780_11585890; do
    if [[ ${GENE} == FBgn0000504 ]]; then
        SYMBOL=dsx
        STRAND=-
        START=""
        END=""
    elif [[ ${GENE} == FBgn0005633 ]]; then
        SYMBOL=fln
        STRAND=-
        START=19989225
        END=19991545
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0000047 ]]; then
        SYMBOL=Act88F
        STRAND=+
        START=15439320
        END=15442840
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0004652 ]]; then
        SYMBOL=fru
        STRAND=-
#        START=
#        END=
        START=""
        END=""
    elif [[ ${GENE} == FBgn0005617 ]]; then
        SYMBOL=msl-1
        STRAND=-
        START=18683620
        END=18688920
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0005616 ]]; then
        SYMBOL=msl-2
        STRAND=-
#        START=3461208
#        END=3467128
        START=3462200
        END=3466200
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0002775 ]]; then
        SYMBOL=msl-3
        STRAND=-
        START=7123033
        END=7125833
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0014340 ]]; then
        SYMBOL=mof
        STRAND=-
        START=5873539
        END=5877039
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0002774 ]]; then
        SYMBOL=mle
        STRAND=-
        START=5973934
        END=5981134
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0019661 ]]; then
        SYMBOL=roX1
        STRAND=-
        START=3856730
        END=3863197
#        START=""
#        END=""
    elif [[ ${GENE} == X_11576780_11585890 ]]; then
        SYMBOL=roX2
        STRAND=+
        START=11578565
        END=11580932
#        START=""
#        END=""
    fi
    ## Remove all PacBio annotations
    ## Modify window to specified positions above
    if [[ ${START} == "" ]]; then
        awk -F "," '$5!~"PB:"' ${INPUT}/${GENE}_annotations.csv > ${ROZ}/${GENE}_annotations.csv
        ## Remove all ChIP values and change columns names for legend
        awk -F "," '{if(NR==1){$3="Female"; $4="Male"}print $1","$2","$3","$4}'\
             ${INPUT}/${GENE}_count_data.csv > ${ROZ}/${GENE}_count_data.csv
    else
        awk -F "," -v start=${START} -v end=${END} \
            '$5!~"PB:" {OFS=","; if($5=="genebody_plus_flank"){\
            $3=start;$4=end} if(NR==1||($3>=start && $4<=end)){print $0}}' \
            ${INPUT}/${GENE}_annotations.csv > ${ROZ}/${GENE}_annotations.csv
        ## Remove all ChIP values and change columns names for legend
        ## Also restrict start and end
        awk -F "," -v start=${START} -v end=${END} \
            'NR==1 || ($2>=start && $2<=end) {if(NR==1)\
            {$3="Female"; $4="Male"}print $1","$2","$3","$4}'\
             ${INPUT}/${GENE}_count_data.csv > ${ROZ}/${GENE}_count_data.csv
    fi
    Rscript ${PROJ}/scripts/integrated_chip_rna/wiggleplots_RNAseq_ChIP_female_male_02avn.R \
        ${ROZ}/${GENE}_annotations.csv \
        ${ROZ}/${GENE}_count_data.csv \
        ${GENE}"("${STRAND}")" \
        ${OUTPUT}/${SPECIES}_${SYMBOL}_${GENE}.png
done

## sim wiggle plots first
SPECIES=sim
for GENE in FBgn0190947 FBgn0183970 FBgn0012823 FBgn0043395 FBgn0268295 FBgn0187546 FBgn0184812 FBgn0269473 FBgn0182116 Scf_X_3403220_3464520 Scf_X_10876200_10885400; do
    if [[ ${GENE} == FBgn0190947 ]]; then
        SYMBOL=GD19448_dsx_ortho
        STRAND=-
        START=""
        END=""
    elif [[ ${GENE} == FBgn0183970 ]]; then
        SYMBOL=GD12234_fln_ortho
        STRAND=-
        START=19473595
        END=19475915
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0012823 ]]; then
        SYMBOL=Act88F
        STRAND=-
        START=9960360
        END=9963880
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0043395 ]]; then
        SYMBOL=fru
        STRAND=+
#        START=
#        END=
        START=""
        END=""
    elif [[ ${GENE} == FBgn0268295 ]]; then
        SYMBOL=GD27005_msl-1_ortho
        STRAND=-
        START=18200305
        END=18205605
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0187546 ]]; then
        SYMBOL=GD15895_msl-2_ortho
        STRAND=-
#        START=3320971
#        END=3326891
        START=3322000
        END=3325850
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0184812 ]]; then
        SYMBOL=GD13090_msl-3_ortho
        STRAND=-
        START=6981986
        END=6984786
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0269473 ]]; then
        SYMBOL=GD28183_mof_ortho
        STRAND=-
        START=5398166
        END=5401666
#        START=""
#        END=""
    elif [[ ${GENE} == FBgn0182116 ]]; then
        SYMBOL=GD10344_mle_ortho
        STRAND=-
        START=2814039
        END=2821239
#        START=""
#        END=""
    elif [[ ${GENE} == Scf_X_3403220_3464520 ]]; then
        SYMBOL=roX1
        STRAND=-
#        START=3450129
        START=3448477
        END=3455044
#        START=""
#        END=""
    elif [[ ${GENE} == Scf_X_10876200_10885400 ]]; then
        SYMBOL=roX2
        STRAND=+
        START=10878030
        END=10880355
#        START=""
#        END=""
    fi
    ## Remove all PacBio annotation
    ## Modify window to specified positions above
    if [[ ${START} == "" ]]; then
        awk -F "," '$5!~"PB:"' ${INPUT}/${GENE}_annotations.csv > ${ROZ}/${GENE}_annotations.csv
        ## Remove all ChIP values and change columns names for legend
        awk -F "," '{if(NR==1){$3="Female"; $4="Male"}print $1","$2","$3","$4}'\
             ${INPUT}/${GENE}_count_data.csv > ${ROZ}/${GENE}_count_data.csv
    else
        awk -F "," -v start=${START} -v end=${END} \
            '$5!~"PB:" {OFS=","; if($5=="genebody_plus_flank"){\
            $3=start;$4=end} if(NR==1||($3>=start && $4<=end)){print $0}}' \
            ${INPUT}/${GENE}_annotations.csv > ${ROZ}/${GENE}_annotations.csv
        ## Remove all ChIP values and change columns names for legend
        ## Also restrict start and end
        awk -F "," -v start=${START} -v end=${END} \
            'NR==1 || ($2>=start && $2<=end){if(NR==1)\
            {$3="Female"; $4="Male"}print $1","$2","$3","$4}'\
             ${INPUT}/${GENE}_count_data.csv > ${ROZ}/${GENE}_count_data.csv
    fi
    Rscript ${PROJ}/scripts/integrated_chip_rna/wiggleplots_RNAseq_ChIP_female_male_02avn.R \
        ${ROZ}/${GENE}_annotations.csv \
        ${ROZ}/${GENE}_count_data.csv \
        ${GENE}"("${STRAND}")" \
        ${OUTPUT}/${SPECIES}_${SYMBOL}_${GENE}.png
done

#rm -r ${ROZ}
