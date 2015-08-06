/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.cegsV_by_gene_sbs;
        run;

    /* 
    proc contents data=w_sample;
    run;
    */

/* SEM with Full covariance */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/full_cov_estimates_yp3.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_model_yp3.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp3 
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_estimates_yp3.csv' dbms=csv replace;
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
