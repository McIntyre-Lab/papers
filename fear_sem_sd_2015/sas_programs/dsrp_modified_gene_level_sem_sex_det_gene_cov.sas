/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sex_det_sbs_gene_level_sym;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

    /* 
    proc contents data=w_sample;
    run;
    */

/* SEM Modified Gene Covariance */
    * The modified gene covarance matrix allows snf, spf45 and fl(2)d to covary
    * with tra2.
    ;
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/modified_gene_cov_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_model.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru,
            tra2 tra her ix -> Yp2,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,

            Spf45 <-> vir =0,

            tra2 <-> vir =0

        ;
        ods output PATHListStd=pathstd ;
        ods output RawRes=resid ;
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
        set pathstd2 ;
        run;

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

    proc export data=resid outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_residuals.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean up */
    proc datasets nolist;
        delete combined;
        delete resid;
        delete pathstd;
        delete pathstd2;
        delete w_sample;
        delete fitstat;
        run; quit;
