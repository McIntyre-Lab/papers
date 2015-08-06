/*******************************************************************************
* Filename: validation_grn_expanstion_dsx_table_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create table of the paper that contains the genes that overlap
* between the 754 added genes and the genes from dsxNull.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data valid;
    set SEM.validation_set_BIC12;
    run;

proc sort data=valid;
    by primary_fbgn BIC;
    run;

data valid2;
    set valid;
    by primary_fbgn;
    if first.primary_fbgn;
    drop model modelnum path;
    run;

proc sort data=valid2 nodupkey;
    by primary_fbgn;
    run;

/* dsx GRN table */
data dsx;
    set valid2;
    where flag_expanded = 1 and flag_dsxnullf_repressed = 1;
    rename flag_dalton_fru_bs = fru_binding_site;
    rename flag_luo_binding_site = dsx_binding_site;
    rename flag_chang_tra_bs = tra_binding_site;
    drop Psex Psex_by_probe bias_toward flag_chang_ds_tra
    flag_chang_female_bias flag_chang_male_bias flag_chang_tra_bias
    flag_chang_traf_bias flag_chang_ups_tra flag_damid_sig
    flag_dsxCtrl_female_induced flag_dsxCtrl_female_repressed
    flag_dsxNullf_induced flag_dsxNullf_repressed flag_goldman_ds_tra
    flag_goldman_female_bias flag_goldman_male_bias flag_goldman_not_ds_dsx
    flag_expanded flag_most_likely
    ;
    run;

proc sort data=dsx;
    by BIC;
    run;

proc export data=dsx outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_overlap_w_dsxNull_repressed_w_luo_go_BIC12.csv' dbms=csv replace;
    putnames=yes;
    run;
