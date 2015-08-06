#!/bin/bash

MCLAB=/home/jfear/mclab
POS_FILE=$MCLAB/cegs_sem_sd_paper/documentation/networks/sex_det_node_positions_w_isoforms.csv
SEM_OUT=$MCLAB/cegs_sem_sd_paper/analysis_output/sem
SCRIPT=$MCLAB/cegs_sem_sd_paper/scripts/sem_graph_sex_det_dsrp.py

# Unconstrained
python $SCRIPT --pos $POS_FILE \
    --sem $SEM_OUT/unconstrained_estimates.csv \
    --out $SEM_OUT/unconstrained_estimates_graph.png \
    --log $SEM_OUT/logs/unconstrained_estimates_graph.log 

# Intra gene exogenous isoforms constrained
python $SCRIPT --pos $POS_FILE \
    --sem $SEM_OUT/constrained_estimates.csv \
    --out $SEM_OUT/constrained_estimates_graph.png \
    --log $SEM_OUT/logs/constrained_estimates_graph.log 

# Intra gene exogenous isoforms partially constrained based on significance in
# unconstrained
python $SCRIPT --pos $POS_FILE \
    --sem $SEM_OUT/partially_constrained_estimates.csv \
    --out $SEM_OUT/partially_constrained_estimates_graph.png \
    --log $SEM_OUT/logs/partially_constrained_estimates_graph.log 

# Spf45 isoforms constrained to 0
python $SCRIPT --pos $POS_FILE \
    --sem $SEM_OUT/spf_isoforms_constrained_estimates.csv \
    --out $SEM_OUT/spf_isoforms_constrained_estimates_graph.png \
    --log $SEM_OUT/logs/spf_isoforms_constrained_estimates_graph.log 

# Intra gene exogenous isoforms constrained and Spf45 isoforms constrained to 0
python $SCRIPT --pos $POS_FILE \
    --sem $SEM_OUT/constrained_w_spf_estimates.csv \
    --out $SEM_OUT/partially_constrained_w_spf_estimates_graph.png \
    --log $SEM_OUT/logs/partially_constrained_w_spf_estimates_graph.log 
