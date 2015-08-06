/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_links_simulations.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: List of SAS programs used for the adding gene analysis. 
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* DO STUFF ON HPC WITH PYTHON AND SAS SEE DOCUMENTATION */
    * I have a python script that generates sas code for all of the
    * different models. I then use HPC to run all of the SAS programs and
    * output a sas_data set for each gene.
    * 
    * SEE: CEGS_sem_adding_links_simulation.xls for more information
    *
    * DATASET: addlink.*
    ;

/* Combine BIC from Adding Links */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addlink '!HOME/tmp/cegs_adding_links_simulation/&i./sas_data'
    * 
    * INPUT: addlink.*
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegs_adding_links_simulation_bic_dist.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_links_simulation.sas';
