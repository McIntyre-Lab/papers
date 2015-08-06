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
    libname dmel '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);


/* DO STUFF ON HPC WITH PYTHON AND SAS SEE DOCUMENTATION */
    * I have a python script that generates sas code for all of the
    * different models. I then use HPC to run all of the SAS programs and
    * output a sas_data set for each gene.
    * 
    * SEE: CEGS_sem_adding_genes.xls for more information
    *
    * DATASET: addgen.*
    ;

/* Combine BIC from Adding genes with Yp2 */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addgen '!MCLAB/cegs_sem_sd_paper/adding_genes/cegsV_nofilter_yp2/sas_data'
    * 
    * INPUT: addgen.*
    *        SEM.cegsV_gene_list
    *
    * DATASET: SEM.cegsV_ag_sex_det_genes_stack_bic
    *          SEM.cegsV_ag_sex_det_genes_sbs_bic
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_sex_det_genes_nofilter_yp2.sas';

/* Look at DSX and see where it fits */
    * I am focusing on DSX for the paper, so I need to look at where it fits in
    * best.
    *
    * INPUT: SEM.cegsV_ag_sex_det_genes_stack_bic
    *
    * DATASET:
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_sex_det_explore_dsx.sas';
