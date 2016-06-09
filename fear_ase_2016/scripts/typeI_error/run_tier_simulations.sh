#!/bin/bash - 
#===============================================================================
#   DESCRIPTION: Run the various TIER simulations for r377, r332, and r365
#===============================================================================

set -o nounset                              # Treat unset variables as an error

PROJ=$MCLAB/cegs_ase_paper
#for LINE in r361 r332 r365 r101; do
for LINE in r101; do
    # Run the three basic models with qTrue = 0.45
    QTRUE=45

    # Run No AI and No Bias simulation, qTRUE is always 0.5 here
    oname=${LINE}_noai_nobias
    python tier_data_simulation.py \
        --input ${PROJ}/pipeline_output/typeI_error/input/${LINE}_RNA_sim_DNA_cnts.csv \
        --output ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.csv \
        --fig ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.pdf


    # Run No AI with Bias simulation
    oname=${LINE}_noai_bias_qTrue0${QTRUE}
    python tier_data_simulation.py \
        --input ${PROJ}/pipeline_output/typeI_error/input/${LINE}_RNA_sim_DNA_cnts.csv \
        --output ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.csv \
        --fig ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.pdf \
        --bias \
        --qtrue 0.${QTRUE} \
        --half


    # Run AI with Bias simulation
    oname=${LINE}_ai_bias_qTrue0${QTRUE}
    python tier_data_simulation.py \
        --input ${PROJ}/pipeline_output/typeI_error/input/${LINE}_RNA_sim_DNA_cnts.csv \
        --output ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.csv \
        --fig ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.pdf \
        --bias \
        --AI \
        --qtrue 0.${QTRUE} \
        --half

    # Run the misspecification models with a range of qTrue
    for QT in 35 375 4 425 45 46 47 48 49 50 51 52 53 54 55 575 6 625 65; do
        oname=${LINE}_misspecification_qTrue0${QT}
        python tier_data_simulation.py \
            --input ${PROJ}/pipeline_output/typeI_error/input/${LINE}_RNA_sim_DNA_cnts.csv \
            --output ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.csv \
            --fig ${PROJ}/pipeline_output/typeI_error/output/${oname}_sim.pdf \
            --bias \
            --qtrue 0.${QT}
    done
done
