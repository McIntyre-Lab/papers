/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_genes_downstream_dsx_BIC12.sas
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

/* Identify genes added downstream of DSX */
    * Using the dsxNull dataset as a valdiation. Need list of genes that
    * were added downstream of DSX.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *        SEM.cegsV_gene_list
    *
    * DATASET: SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_flag_ds_dsx_BIC12.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_dsx_gene_list_nofilter_yp2_BIC12.sas';

/* Compare DSX list to Luo 2011 */
    * Luo 2011 identified 58 genes with dsx binding sites in them. I want
    * to compare their list to the genes added downstream of dsx.
    *
    * INPUT: SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12
    * 
    * DATASET: SEM.luo2011
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_dsx_compare_luo2011_BIC12.sas';

/* Create list of genes that fit downstream of DSX the best */
    * I am wanting to see what happens when I try to add multiple genes
    * simultaneously downstream of DSX, so I need to create a list of genes
    * that had their best fit downstream of DSX (model 3: Dsx -> gene).
    *
    * INPUT: SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12
    * 
    * DATASET:
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_dsx_add_multiple_downstream_BIC12.sas';
