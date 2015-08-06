/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sex_det_sbs_combine_sym;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

    /* 
    proc contents data=w_sample;
    run;
    */

/* SEM with Full covariance */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/fullcov_isoform_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/fullcov_isoform_model.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
        run;

    proc printto LOG=LOG PRINT=PRINT;
    run;

    data pathstd2;
        set pathstd;
        drop arrow;
        run;

    data combined;
        retain SpecType Var1 Var2 Parameter Estimate StdErr tValue;
        length Var2 $8.;
        set pathstd2 covstd;
        run;

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/fullcov_isoform_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean up */
    proc datasets nolist;
        delete combined;
        delete covstd;
        delete pathstd;
        delete pathstd2;
        delete w_sample;
        run; quit;
