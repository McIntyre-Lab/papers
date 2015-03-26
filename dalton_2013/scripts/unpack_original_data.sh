#!/bin/bash
#
#PBS -q submit
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#

work=/scratch/hpc/mccrory/fru_network

cd ${work}/original_data/
tar -zxvf Data_${date}.tgz
