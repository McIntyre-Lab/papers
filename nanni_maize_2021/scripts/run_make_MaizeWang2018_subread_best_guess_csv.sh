#!/bin/sh

## Set paths
PROJ=/blue/mcintyre/share/maize_ainsworth
DATA=/blue/mcintyre/share/transcript_distance/MaizeWang2018
DESIGN=$PROJ/design_files
OUTFILE=${DESIGN}/MaizeWang2018_subread_best_guess.csv

echo "fullpath,filename,UFbestGuess_LID,UFbestGuess_pool,UFbestGuess_insertSize,UFbestGuess_barcodeID,UFbestGuess_species,UFbestGuess_tissue" > ${OUTFILE}


for FULLPATH in $(find $DATA -name '*.subreads.fastq'); do
    FILE=$(basename ${FULLPATH})
    LIDGUESS=none
    ## Check for pool 1 LID values from excel file provided by Bo Wang
    ## Some LID values had a "2" in the upper directory levels but not in the excel
    ##     or in the SMRT50 directory (both are checked)
    for LID in LID50283 LID50284 LID50285 LID502302 LID502303 LID50302 LID50303; do
        ## Check if LID is in the full path surrounded by "/"
        if [[ "${FULLPATH}" =~ "/${LID}/" ]]; then
            LIDGUESS=${LID}
            POOLGUESS=1
            BCGUESS="BC1|BC2|BC3|BC4|BC5|BC6|BC7|BC8|BC9"
            SPGUESS="maize|sorghum"
            TGUESS="maize_silk|maize_20_DAP_pericarp|maize_bract|maize_seedling_shoot|maize_14_DAG_seedling|maize_14_DAG_seedling_leaf|sorgum_14_DAG_seedling|sorgum_seedling_shoot|sorgum_20_DAG_endosperm"
        fi
    done
    ## Check for pool 2 LID values from excel file provided by Bo Wang
    ## Some LID	values had a "2" in the	upper directory	levels but not in the excel
    ## 	   or in the SMRT50 directory (both are	checked)
    for LID in LID50288 LID50289 LID50290 LID502304 LID502305 LID50304 LID50305; do
        ## Check if LID is in the full path surrounded by "/"
        if [[ "${FULLPATH}" =~ "/${LID}/" ]]; then
            LIDGUESS=${LID}
            POOLGUESS=2
       	    BCGUESS="BC1|BC2|BC3|BC4|BC5|BC6|BC7|BC8|BC9"
            SPGUESS="sorghum"
            TGUESS="sorgum_14_DAG_seedling_root|sorgum_14_DAG_seedling_leaf|sorgum_integument|sorgum_1-2_cm_inflorescence|sorgum_0.5-1_cm_inflorescence|sorgum_pollen|sorgum_.5_cm_inflorescence|sorgum_20_DAP_embryo|sorgum_Jura_7"
        fi
    done
    ## Set insert sizes for pools from main LID directories
    if [[ ${LIDGUESS} != "none" ]]; then
        if [[ ${LIDGUESS} == "LID50283" || ${LIDGUESS} == "LID50288" ]]; then
            INSGUESS="1kb"
        elif [[ ${LIDGUESS} == "LID50284" || ${LIDGUESS} == "LID50289" ]]; then
            INSGUESS="1to2kb"
        elif [[ ${LIDGUESS} == "LID50285" || ${LIDGUESS} == "LID50290" ]]; then
       	    INSGUESS="2to3kb"
        elif [[ ${LIDGUESS} == "LID502302" || ${LIDGUESS} == "LID50302" || ${LIDGUESS} == "LID502304" || ${LIDGUESS} == "LID50304" ]]; then
            INSGUESS="3to5kb"
        elif [[ ${LIDGUESS} == "LID502303" || ${LIDGUESS} == "LID50303" || ${LIDGUESS} == "LID502305" || ${LIDGUESS} == "LID50305" ]]; then
            INSGUESS="over5kb"
        else
            echo "ERROR LABELLING ${FULLPATH}"
            break
        fi
    ## Set values for the "maize" directory of insert sizes
    else
        if [[ "${FULLPATH}" =~ "/maize/" ]]; then
            POOLGUESS=1
       	    BCGUESS="BC1|BC2|BC3|BC4|BC5|BC6"
            SPGUESS="maize"
       	    TGUESS="maize_silk|maize_20_DAP_pericarp|maize_bract|maize_seedling_shoot|maize_14_DAG_seedling|maize_14_DAG_seedling_leaf"
            ## Get insert size from name
            if [[ "${FULLPATH}" =~ "/BC_under1kb/" ]]; then
       	       	INSGUESS="under1kb"
            elif [[ "${FULLPATH}" =~ "/BC_1to2kb/" ]]; then
                INSGUESS="1to2kb"
       	    elif [[ "${FULLPATH}" =~ "/BC_2to3kb/" ]]; then
       	       	INSGUESS="2to3kb"
       	    elif [[ "${FULLPATH}" =~ "/BC_3to5kb/" ]]; then
       	       	INSGUESS="3to5kb"
       	    elif [[ "${FULLPATH}" =~ "/BC_4to6kb/" ]]; then
       	       	INSGUESS="4to6kb"
       	    elif [[ "${FULLPATH}" =~ "/BC_5to10kb/" ]]; then
       	       	INSGUESS="5to10kb"
            else
                echo "ERROR LABELLING ${FULLPATH}"
                break
            fi
        else
            echo "ERROR LABELLING ${FULLPATH}"
            break
        fi
    fi
    echo "${FULLPATH},${FILE},${LIDGUESS},${POOLGUESS},${INSGUESS},${BCGUESS},${SPGUESS},${TGUESS}" >> ${OUTFILE}
done
