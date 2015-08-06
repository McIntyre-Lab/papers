/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* How many fusions are there for each gene? */
    data sex;
        set SEM.cegsV_sex_det_stack;
        keep symbol_cat fusion_id ;
        run;

    proc sort data=sex nodupkey;
        by symbol_cat fusion_id;
        run;

    proc freq data=sex;
        table symbol_cat;
        run;

/* Are the fusions normally distributed? */
    proc sort data=SEM.cegsV_sex_det_stack;
        by fusion_id;
        run;

    proc univariate data=SEM.cegsV_sex_det_stack normal plots;
        by fusion_id;
        var mean_exp;
        ods output TestsForNormality=normtest;
        run;

    data flag_fail_normality;
        set normtest;
        where test="Shapiro-Wilk";
        if pValue = . then flag_fail_normality = .;
        else if pValue <= 0.05 then flag_fail_normality = 1;
        else flag_fail_normality = 0;
        keep fusion_id flag_fail_normality;
        run;

    proc freq data=flag_fail_normality;
        table flag_fail_normality;
        run; * 33 pass, 67 fail;

/* Do I lose any genes if I drop fusions that fail normality? */
    data pass;
        set flag_fail_normality;
        where flag_fail_normality = 0;
        run;

    proc sort data=pass;
        by fusion_id;
        run;

    proc sort data=sex;
        by fusion_id;
        run;

    data merged;
        merge sex (in=in1) pass (in=in2);
        by fusion_id;
        if in2;
        run;

    proc freq data=merged;
        table symbol_cat;
        run; * lost fl_2_d, snf, tra2;

/* Look at tra2 closer */
    proc univariate data=SEM.cegsV_sex_det_stack(where=(symbol_cat eq 'tra2')) normal plots;
        by fusion_id;
        var mean_exp;
        run; * F40928_SI is the most normal fusion.;

/* Look at snf closer */
    proc univariate data=SEM.cegsV_sex_det_stack(where=(symbol_cat eq 'snf')) normal plots;
        by fusion_id;
        var mean_exp;
        run; * S13633_SI is the most normal fusion.;

/* Look at fl_2_d closer */
    proc univariate data=SEM.cegsV_sex_det_stack(where=(symbol_cat eq 'fl_2_d')) normal plots;
        by fusion_id;
        var mean_exp;
        run; * F40441_SI is the most normal fusion.;

/* Clean Up */
proc datasets;
    delete flag_fail_normality;
    delete merged;
    delete normtest;
    delete pass;
    delete sex;
    run;
