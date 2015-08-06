/*******************************************************************************
* Filename: enrichment_list_enrichment.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Using a variety of gene lists test for enrichment.
*
*******************************************************************************/

/* Libraries
    libname SEM '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
    libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/*
flag_expanded
flag_most_likely
flag_dsx_exp
flag_dsx_exp_most_likely
flag_fru_exp
flag_fru_exp_most_likely
modelnum
path
flag_dalton_fru_bs
flag_dsxCtrl_female_induced
flag_dsxCtrl_female_repressed
flag_dsxNullf_induced
flag_dsxNullf_repressed
flag_damid_sig
flag_luo_binding_site
flag_chang_male_bias
flag_chang_female_bias
flag_chang_tra_bias
flag_chang_traf_bias
flag_chang_ds_tra
flag_chang_ups_tra
flag_chang_tra_bs
bias_toward
Psex
Psex_by_probe
flag_mcintyre_sex_bias_splice
flag_goldman_female_bias
flag_goldman_male_bias
flag_goldman_ds_tra
flag_goldman_not_ds_dsx
*/

proc sort data=SEM.validation_set_bic12;
    by primary_fbgn BIC;
    run;

data collapse;
    set SEM.validation_set_bic12;
    by primary_fbgn;
    if first.primary_fbgn;
    if flag_chang_ds_tra = '.' then flag_chang_ds_tra  = 0;
    if flag_chang_female_bias = '.' then flag_chang_female_bias  = 0;
    if flag_chang_male_bias = '.' then flag_chang_male_bias  = 0;
    if flag_chang_tra_bias = '.' then flag_chang_tra_bias  = 0;
    if flag_chang_tra_bs = '.' then flag_chang_tra_bs  = 0;
    if flag_chang_traf_bias = '.' then flag_chang_traf_bias  = 0;
    if flag_chang_ups_tra = '.' then flag_chang_ups_tra  = 0;
    if flag_dalton_fru_bs = '.' then flag_dalton_fru_bs  = 0;
    if flag_damid_sig = '.' then flag_damid_sig  = 0;
    if flag_dsxCtrl_female_induced = '.' then flag_dsxCtrl_female_induced  = 0;
    if flag_dsxCtrl_female_repressed = '.' then flag_dsxCtrl_female_repressed  = 0;
    if flag_dsxNullf_induced = '.' then flag_dsxNullf_induced  = 0;
    if flag_dsxNullf_repressed = '.' then flag_dsxNullf_repressed  = 0;
    if flag_goldman_ds_tra = '.' then flag_goldman_ds_tra  = 0;
    if flag_goldman_female_bias = '.' then flag_goldman_female_bias  = 0;
    if flag_goldman_male_bias = '.' then flag_goldman_male_bias  = 0;
    if flag_goldman_not_ds_dsx = '.' then flag_goldman_not_ds_dsx  = 0;
    if flag_luo_binding_site = '.' then flag_luo_binding_site  = 0;
    if flag_mcintyre_sex_bias_splice = '.' then flag_mcintyre_sex_bias_splice  = 0;
    if flag_dsxNullf_repressed = 1 and flag_expanded = 1 then flag_dsx_overlap = 1;
    else flag_dsx_overlap = 0;
    drop model BIC path modelnum;
    run;

/* Sanity check */
proc freq data=collapse;
    tables flag_expanded;
    run; * 754 ok;

proc freq data=collapse;
    tables flag_dsx_overlap;
    run; * 754 ok;

/* Freqs */
    proc freq data=collapse;
        tables flag_expanded*flag_chang_ds_tra /chisq exact expect;
        run; *no enrich;

    proc freq data=collapse;
        tables flag_expanded*flag_dsxNullf_repressed /chisq exact expect;
        run; *no enrich;

    proc freq data=collapse;
        tables flag_expanded*flag_chang_tra_bs /chisq exact expect;
        run; *no enrich;

    proc freq data=collapse;
        tables flag_expanded*flag_mcintyre_sex_bias_splice /chisq exact expect;
        run; *enrich;

    /* DSX */
        proc freq data=collapse;
            tables flag_dsx_exp*flag_luo_binding_site /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp_most_likely*flag_luo_binding_site /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_overlap*flag_luo_binding_site /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp*flag_dsxNullf_induced /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp_most_likely*flag_dsxNullf_induced /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp*flag_dsxNullf_repressed /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp_most_likely*flag_dsxNullf_repressed /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp*flag_chang_ds_tra /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp_most_likely*flag_chang_ds_tra /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp*flag_goldman_ds_tra /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsx_exp_most_likely*flag_goldman_ds_tra /chisq exact expect;
            run; *no enrich;

        proc freq data=collapse;
            tables flag_dsxCtrl_female_induced*flag_chang_female_bias /chisq exact expect;
            run; *enrich;

        proc freq data=collapse;
            tables flag_dsxCtrl_female_induced*flag_goldman_female_bias /chisq exact expect;
            run; *enrich;


    /* FRU */
        proc freq data=collapse;
            tables flag_fru_exp*flag_dalton_fru_bs /chisq exact expect;
            run; *no enrich but close;
    
        proc freq data=collapse;
            tables flag_fru_exp_most_likely*flag_dalton_fru_bs /chisq exact expect;
            run; *no enrich;


/* Misc stuff */
data tmp;
    set collapse;
    if flag_chang_ds_tra = 1 or flag_chang_ups_tra = 1 then flag_tra =1;
    else flag_tra = 0;
    run;

proc freq data=tmp;
    tables flag_expanded*flag_tra;
    run;



