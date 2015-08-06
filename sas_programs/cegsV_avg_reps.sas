/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname norm '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab virgins from CEGS normalized and centered data */
data virgin;
    set NORM.ccfus_norm_centered;
    where mating_status eq 'V';
    keep fusion_id line rep uq_log_uq_center;
    run;

proc sort data=virgin;
    by fusion_id line ;
    run;

/* Average across replicates */
proc means data=virgin noprint;
    by fusion_id line ;
    output out=means mean(uq_log_uq_center)=mean_exp;
    run;

/* Grab annotations from DMEL */
data anno;
    set DMEL551.fb551_si_fusions_unique_flagged;
    keep fusion_id genes_per_fusion symbol_cat FBgn_cat;
    run;

proc sort data=anno;
    by fusion_id;
    run;

proc sort data=means;
    by fusion_id;
    run;

data merged;
    merge means (in=in1) anno (in=in2);
    by fusion_id;
    if in1;
    drop _type_ _freq_;
    run;

data SEM.cegs_virgin_norm_cent;
    set merged;
    run;

/* Clean Up */
proc datasets nolist;
    delete virgin;
    delete means;
    delete anno;
    delete merged;
    run; quit;
