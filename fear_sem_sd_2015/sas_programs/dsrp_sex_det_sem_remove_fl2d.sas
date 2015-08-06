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

/* BASELINE Gene Covaraince */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_fl2d.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_fl2d.lst' NEW;
    run;

    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            Spf45 snf vir -> Sxl_A,
            vir Sxl_A -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp2,

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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_fl2d.csv' dbms=csv replace;
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
