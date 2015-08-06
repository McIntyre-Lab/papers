/*******************************************************************************
* Filename: misc_test_gene_expression.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Figure out how InR and the sex determination gene expression
* compares with all of the other genes.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data for_plot;
    set SEM.cegsv_by_gene_stack;
    if sym eq 'Sxl' or
       sym eq 'fl_2_d' or 
       sym eq 'dsx' or 
       sym eq 'fru' or 
       sym eq 'tra' or 
       sym eq 'tra2' or 
       sym eq 'Spf45' or 
       sym eq 'snf' or 
       sym eq 'ix' or 
       sym eq 'her' or 
       sym eq 'Yp1' or 
       sym eq 'Yp2' or 
       sym eq 'Yp3' or 
       sym eq 'vir' then flag_sex = 1;
    else if sym eq 'FBgn0013984' then flag_sex = 2;
    else flag_sex = 0;
    run;

libname tmp '/home/jfear/tmp';
data tmp.for_plot;
    set for_plot;
    run;
