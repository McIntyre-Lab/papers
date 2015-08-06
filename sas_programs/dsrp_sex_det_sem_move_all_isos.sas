/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
    libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sas_data';
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

/* BASELINE SEM with constrained exogenous GENES */
    /*
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2,

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
    */

/* Remove Spf45 from model */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_spf45_estimates_full.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_spf45_estimates_full.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
        run;

    * print log and output back to listings;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_spf45_estimates_full.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Move Spf45 To Bottom Model */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_spf45_bottom_estimates_full.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_bottom_estimates_full.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            Spf45 tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

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

    * print log and output back to listings;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_bottom_estimates_full.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Add Spf45 link to fru, suggested by residuals */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/spf45_add_fru_estimates_full.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/spf45_add_fru_estimates_full.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra Spf45 -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

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

    * print log and output back to listings;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/spf45_add_fru_estimates_full.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Move vir To Tra2 Model */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_vir_to_tra2_estimates_full.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_estimates_full.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B Sxl_A Sxl_B Sxl_C -> tra,
            vir -> Sxl_A tra2_A tra2_B,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3 ,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> vir =0,

            Spf45 <-> vir =0

        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
        run;

    * print log and output back to listings;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_estimates_full.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Move Yp2 to top Model */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_yp2_to_top_estimates_full.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_yp2_to_top_estimates_full.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            Yp1 Yp2 Yp3 -> fl_2_d_A fl_2_d_B,
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,

            Yp1 <-> her  =0,
            Yp1 <-> ix  =0,
            Yp1 <-> snf  =0,
            Yp1 <-> Spf45 =0,
            Yp1 <-> tra2_A  =0,
            Yp1 <-> tra2_B  =0,
            Yp1 <-> vir =0,

            Yp2 <-> her  =0,
            Yp2 <-> ix  =0,
            Yp2 <-> snf  =0,
            Yp2 <-> Spf45 =0,
            Yp2 <-> tra2_A  =0,
            Yp2 <-> tra2_B  =0,
            Yp2 <-> vir =0,

            Yp3 <-> her  =0,
            Yp3 <-> ix  =0,
            Yp3 <-> snf  =0,
            Yp3 <-> Spf45 =0,
            Yp3 <-> tra2_A  =0,
            Yp3 <-> tra2_B  =0,
            Yp3 <-> vir =0,

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

    * print log and output back to listings;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_yp2_to_top_estimates_full.csv' dbms=csv replace;
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
