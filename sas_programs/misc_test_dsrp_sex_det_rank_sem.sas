/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sexdet_comb_sym_from_rank;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* (1) SEM with unconstrained covariances*/
    options ps = 70 ls=80;
    proc calis data=w_sample method=ml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45_A Spf45_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/rank_sem_output/unconstrained_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

/* (2) SEM with covarainces between exogenous GENES constrained to 0 */
    options ps = 70 ls=80;
    proc calis data=w_sample method=ml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45_A Spf45_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45_A =0,
            fl_2_d_A <-> Spf45_B  =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45_A =0,
            fl_2_d_B <-> Spf45_B  =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45_A =0,
            her <-> Spf45_B  =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45_A =0,
            ix <-> Spf45_B  =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45_A =0,
            snf <-> Spf45_B  =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45_A <-> tra2_A  =0,
            Spf45_A <-> tra2_B  =0,
            Spf45_A <-> vir =0,

            Spf45_B <-> tra2_A  =0,
            Spf45_B <-> tra2_B  =0,
            Spf45_B <-> vir =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/rank_sem_output/constrained_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

/* (3) SEM with covarainces between select GENES constrained to 0 */
    options ps = 70 ls=80;
    proc calis data=w_sample method=ml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45_A Spf45_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,
            fl_2_d_B <--> Spf45_B = 0,
            her <--> Spf45_A = 0,
            her <--> Spf45_B = 0,
            ix <--> fl_2_d_B = 0,
            snf <--> Spf45_A = 0,
            snf <--> Spf45_B = 0,
            snf <--> fl_2_d_B = 0,
            snf <--> her = 0,
            snf <--> ix = 0,
            tra2_A <--> Spf45_A = 0,
            tra2_A <--> Spf45_B = 0,
            tra2_A <--> fl_2_d_A = 0,
            tra2_A <--> fl_2_d_B = 0,
            tra2_B <--> Spf45_A = 0,
            tra2_B <--> Spf45_B = 0,
            tra2_B <--> fl_2_d_A = 0,
            tra2_B <--> fl_2_d_B = 0,
            tra2_B <--> her = 0,
            tra2_B <--> ix = 0,
            tra2_B <--> snf = 0,
            vir <--> Spf45_A = 0,
            vir <--> Spf45_B = 0,
            vir <--> fl_2_d_B = 0,
            vir <--> ix = 0,
            vir <--> snf = 0,
            vir <--> tra2_A = 0,
            vir <--> tra2_B = 0
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/rank_sem_output/partially_constrained_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

/* (4) SEM with covarainces between Spf isoforms constrained to 0 */
    options ps = 70 ls=80;
    proc calis data=w_sample method=ml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45_A Spf45_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

            Spf45_A <-> Spf45_B =0
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/rank_sem_output/spf_isoforms_constrained_estimates.csv' dbms=csv replace;
        putnames=yes;
        run;

/* (5) SEM with covarainces between exogenous GENES constrained to 0 and Spf within gene constrained to 0 */
    options ps = 70 ls=80;
    proc calis data=w_sample method=ml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        path
            fl_2_d_A fl_2_d_B Spf45_A Spf45_B snf vir -> Sxl_A Sxl_B Sxl_C,
            fl_2_d_A fl_2_d_B vir Sxl_A Sxl_B Sxl_C -> tra,
            tra2_A tra2_B tra -> fru,
            tra2_A tra2_B tra her ix -> Yp1 Yp2 Yp3,

            fl_2_d_A <-> her  =0,
            fl_2_d_A <-> ix  =0,
            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45_A =0,
            fl_2_d_A <-> Spf45_B  =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> her  =0,
            fl_2_d_B <-> ix  =0,
            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45_A =0,
            fl_2_d_B <-> Spf45_B  =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

            her <-> ix  =0,
            her <-> snf  =0,
            her <-> Spf45_A =0,
            her <-> Spf45_B  =0,
            her <-> tra2_A  =0,
            her <-> tra2_B  =0,
            her <-> vir =0,

            ix <-> snf  =0,
            ix <-> Spf45_A =0,
            ix <-> Spf45_B  =0,
            ix <-> tra2_A  =0,
            ix <-> tra2_B  =0,
            ix <-> vir =0,

            snf <-> Spf45_A =0,
            snf <-> Spf45_B  =0,
            snf <-> tra2_A  =0,
            snf <-> tra2_B  =0,
            snf <-> vir =0,

            Spf45_A <-> tra2_A  =0,
            Spf45_A <-> tra2_B  =0,
            Spf45_A <-> vir =0,

            Spf45_B <-> tra2_A  =0,
            Spf45_B <-> tra2_B  =0,
            Spf45_B <-> vir =0,

            Spf45_A <-> Spf45_B =0,

            tra2_A <-> vir =0,

            tra2_B <-> vir =0
        ;
        ods output PATHListStd=pathstd ;
        ods output PATHCovVarsStd=covstd ;
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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/rank_sem_output/constrained_w_spf_estimates.csv' dbms=csv replace;
        putnames=yes;
        run;

    * Clean up;
    proc datasets nolist;
        delete combined;
        delete covstd;
        delete pathstd;
        delete pathstd2;
        delete w_sample;
        run; quit;
