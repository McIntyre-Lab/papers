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

/* SEM constrained model with all Sxl isoforms */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_all_sxl.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_all_sxl.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C-> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_all_sxl.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_A,B */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlAB.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlAB.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A Sxl_B,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlAB.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_A,C */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlAC.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlAC.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlAC.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_B,C */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlBC.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlBC.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlBC.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_B,C */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlB.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlB.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_B,
            fl_2_d_A fl_2_d_B vir Sxl_B -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlB.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_B */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlB.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlB.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_B,
            fl_2_d_A fl_2_d_B vir Sxl_B -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlB.csv' dbms=csv replace;
    putnames=yes;
    run;

/* SEM constrained model with Sxl_C */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_estimates_sxlC.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlC.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45 <-> tra2_A  =0,
            Spf45 <-> tra2_B  =0,
            Spf45 <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlC.csv' dbms=csv replace;
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
