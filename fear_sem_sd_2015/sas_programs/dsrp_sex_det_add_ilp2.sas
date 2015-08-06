/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Add ILP2 to SEM dataset */
    data ilp2;
        set SEM.dsrp;
        rename CG8167 = ilp2;
        keep patRIL matRIL CG8167;
        run;

    proc sort data=ilp2;
        by patRIL matRIL;
        run;

    proc sort data=SEM.dsrp_sex_det_sbs_combine_sym;
        by patRIL matRIL;
        run;

    data merged;
        merge SEM.dsrp_sex_det_sbs_combine_sym (in=in1) ilp2 (in=in2);
        by patRIL matRIL;
        if in1 and in2;
        run;

/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set merged;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* BASELINE Gene Covaraince */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/add_ilp2_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            ilp2 -> tra2_A tra2_b,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> ilp2  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> ilp2  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> ilp2  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45 =0,
            ix <-> ilp2  =0,
            ix <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> ilp2  =0,
            snf <-> vir =0,

            Spf45 <-> ilp2  =0,
            Spf45 <-> vir =0,

            ilp2 <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;







/* BASELINE Gene Covaraince */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/add_ilp2_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            ilp2 -> tra2_A tra2_b ix,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> ilp2  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> ilp2  =0,
            fl_2_d_B <-> vir =0,

            her <-> snf  =0,
            her <-> Spf45 =0,
            her <-> ilp2  =0,
            her <-> vir =0,

            snf <-> Spf45 =0,
            snf <-> ilp2  =0,
            snf <-> vir =0,

            Spf45 <-> ilp2  =0,
            Spf45 <-> vir =0,

            ilp2 <-> vir =0
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;



/* BASELINE Gene Covaraince */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/add_ilp2_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2 ipl2,

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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/add_ilp2_estimates.csv' dbms=csv replace;
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
