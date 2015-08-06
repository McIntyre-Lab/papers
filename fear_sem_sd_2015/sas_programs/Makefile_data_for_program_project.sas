/*******************************************************************************
* Filename: Makefile_data_for_program_project.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Lauren is working on a large program project grant. These files
* were used to generate some data for that grant.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Summarize added genes upstream of fru and dsx but not DS of SXL */
    * Sergey asked for a list of genes that are putative upstream of fru and
    * dsx, but not downstream of SXL. 
    *
    * I am particularly interested in:
    *   Model_15 (gene -> dsx)
    *   Model_16 (gene -> fru)
    *   Model_26 (tra -> gene -> dsx)
    *   Model_27 (tra -> gene -> fru)
    *   Model_34 (tra2 -> gene -> dsx)
    *   Model_35 (tra2 -> gene -> fru)
    *
    * But not:
    *   Model_1 (Sxl -> gene)
    *
    * INPUT: SEM.cegsv_ag_w_bic_cutoff
    *        dmel.symbol2fbgn
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/program_project_genes_us_fru_and_dsx_not_ds_sxl.csv
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/program_project_upstream_fru_dsx_not_ds_sxl.sas';
