/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_genes_downstream_fru.sas
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

/* Identify genes added downstream of FRU */
    * It would be useful to know how many genes are downstream of fru. 
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *
    * DATASET: SEM.cegsV_ag_yp2_flag_ds_fru
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_flag_ds_fru.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_fru_gene_list_nofilter_yp2.sas';

/* Fru BS analysis */
    * Using the previous Dalton 2013 Fru BS data, are the 9 genes that are best
    * downstream of fru have the fru BS?
    *
    * libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
    *
    * INPUT: SEM.cegsV_ag_yp2_flag_ds_fru
    *        FRU.motif_flags_and_cnts
    *
    * DATASET:
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_fru_bs_analysis.sas';
