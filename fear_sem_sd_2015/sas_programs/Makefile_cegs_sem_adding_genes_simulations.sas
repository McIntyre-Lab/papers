/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_genes.sas
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
    * SEE: CEGS_sem_simulation.xls for more information
    *
    * DATASET: addgen.*
    ;

/* 80 gene simulation */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addgen !HOME/tmp/cegs_adding_genes_simulation/&numgenes._gene_sim/&i./sas_data
    * 
    * INPUT: addgen.*
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegs_adding_genes_simulation_&numgenes._gene_sim_bic_dist.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_80gene_simulation.sas';

/* 800 gene simulation */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addgen !HOME/tmp/cegs_adding_genes_simulation/&numgenes._gene_sim/&i./sas_data
    * 
    * INPUT: addgen.*
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegs_adding_genes_simulation_&numgenes._gene_sim_bic_dist.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_800gene_simulation.sas';

/* 8000 gene simulation */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addgen !HOME/tmp/cegs_adding_genes_simulation/&numgenes._gene_sim/&i./sas_data
    * 
    * INPUT: addgen.*
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegs_adding_genes_simulation_&numgenes._gene_sim_bic_dist.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_800gene_simulation.sas';
