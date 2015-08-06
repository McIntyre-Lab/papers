/*******************************************************************************
* Filename: validation_freqs_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Focusing on DSX models, what are the FREQS for BS and all things
* dsx.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Set missing flags to 0 for chisq tests */
data valid;
    set SEM.validation_set_BIC12;
    run;

data valid2;
    set valid;
    if flag_chang_ds_tra = '.' then flag_chang_ds_tra  = 0;
    if flag_chang_female_bias = '.' then flag_chang_female_bias  = 0;
    if flag_chang_male_bias = '.' then flag_chang_male_bias  = 0;
    if flag_chang_tra_bias = '.' then flag_chang_tra_bias  = 0;
    if flag_chang_tra_bs = '.' then flag_chang_tra_bs  = 0;
    if flag_chang_traf_bias = '.' then flag_chang_traf_bias  = 0;
    if flag_chang_ups_tra = '.' then flag_chang_ups_tra  = 0;
    if flag_luo_binding_site = '.' then flag_luo_binding_site  = 0;
    if flag_damid_sig = '.' then flag_damid_sig  = 0;
    if flag_dsxCtrl_female_induced = '.' then flag_dsxCtrl_female_induced  = 0;
    if flag_dsxCtrl_female_repressed = '.' then flag_dsxCtrl_female_repressed  = 0;
    if flag_dsxNullf_induced = '.' then flag_dsxNullf_induced  = 0;
    if flag_dsxNullf_repressed = '.' then flag_dsxNullf_repressed  = 0;
    if flag_goldman_ds_tra = '.' then flag_goldman_ds_tra  = 0;
    if flag_goldman_female_bias = '.' then flag_goldman_female_bias  = 0;
    if flag_goldman_male_bias = '.' then flag_goldman_male_bias  = 0;
    if flag_goldman_not_ds_dsx = '.' then flag_goldman_not_ds_dsx  = 0;
    if flag_grn_expansion = '.' then flag_grn_expansion  = 0;
    if flag_mcintyre_sex_bias_splice = '.' then flag_mcintyre_sex_bias_splice  = 0;
    drop model bic modelnum path;
    run;
    
/* Merge on dsx SEM model flags  */
proc sort data=valid2;
    by symbol;
    run;

proc sort data=SEM.cegsV_ag_yp2_flag_ds_dsx_bic12;
    by symbol;
    run;

data valid3;
    merge valid2 (in=in1) SEM.cegsV_ag_yp2_flag_ds_dsx_bic12 (in=in2);
    by symbol;
    if in1;
    run;

/* Check genes mostlikely ds dsx */




proc sort data=valid2 nodupkey;
    by primary_fbgn;
    run;

proc sort data=valid2 nodupkey;
    by primary_fbgn;
    run;

proc freq data=valid2;
    tables flag_grn_expansion;
    tables flag_grn_expansion*(flag_dsxNullf_repressed
    flag_chang_ds_tra
    flag_chang_ups_tra
    flag_chang_tra_bs
    flag_goldman_ds_tra) / chisq exact;
    run;

data check;
    set valid2;
    if flag_chang_tra_bs=1 and flag_grn_expansion=1;
    run;

proc freq data=check;
    tables flag_grn_expansion*(flag_dsxNullf_repressed
    flag_chang_ds_tra
    flag_chang_ups_tra
    flag_chang_tra_bs
    flag_goldman_ds_tra);
    run;

