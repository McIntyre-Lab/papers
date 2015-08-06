/********************************************************************************
* I am testing the different ways to specify the model using proc calis.
********************************************************************************/
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Create Input Dataset, combine parental ids for simplicity */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sex_det_sbs_combine_sym;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* PATH specification */
    options ps = 70 ls=80;
    proc printto PRINT='/home/jfear/Desktop/path.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
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
        run;

    proc printto PRINT=PRINT;
    run;

/* Lineqs */
    options ps = 70 ls=80;
    proc printto PRINT='/home/jfear/Desktop/lineqs.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
        lineqs
            Sxl_A = P1 fl_2_d_A + P2 fl_2_d_B + P3 Spf45 + P4 snf + P5 vir + E1,
            tra = P6 fl_2_d_A + P7 fl_2_d_B + P8 vir + P9 Sxl_A + E2,
            fru = P10 tra2_A + P11 tra2_B + P12 tra + E3,
            Yp2 = P13 tra2_A + P14 tra2_B + P15 tra + P16 her + P17 ix + E4;
        cov
            fl_2_d_A her  =0,
            fl_2_d_A ix  =0,
            fl_2_d_A snf  =0,
            fl_2_d_A Spf45 =0,
            fl_2_d_A tra2_A  =0,
            fl_2_d_A tra2_B  =0,
            fl_2_d_A vir =0,
            fl_2_d_B her  =0,
            fl_2_d_B ix  =0,
            fl_2_d_B snf  =0,
            fl_2_d_B Spf45 =0,
            fl_2_d_B tra2_A  =0,
            fl_2_d_B tra2_B  =0,
            fl_2_d_B vir =0,
            her ix  =0,
            her snf  =0,
            her Spf45 =0,
            her tra2_A  =0,
            her tra2_B  =0,
            her vir =0,
            ix snf  =0,
            ix Spf45 =0,
            ix tra2_A  =0,
            ix tra2_B  =0,
            ix vir =0,
            snf Spf45 =0,
            snf tra2_A  =0,
            snf tra2_B  =0,
            snf vir =0,
            Spf45 tra2_A  =0,
            Spf45 tra2_B  =0,
            Spf45 vir =0,
            tra2_A vir =0,
            tra2_B vir =0
        ;
        run;
    proc printto PRINT=PRINT;
    run;

/* lismod */
    options ps = 70 ls=80;
    proc printto PRINT='/home/jfear/Desktop/lismod.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=test;
        lismod 
            yvar = Sxl_A tra fru Yp2,
            xvar = fl_2_d_A fl_2_D_B Spf45 snf vir tra2_A tra2_B her ix;
        matrix _BETA_ [2,1] = beta1,
                      [3,2] = beta2,
                      [4,2] = beta3
        ;
        matrix _GAMMA_ [1,1] = gamma1 ,
                       [1,2] = gamma2 ,
                       [1,3] = gamma3 ,
                       [1,4] = gamma4 ,
                       [1,5] = gamma5 ,
                       [2,1] = gamma6 ,
                       [2,2] = gamma7 ,
                       [2,5] = gamma8 ,
                       [3,6] = gamma9 ,
                       [3,7] = gamma10,
                       [4,6] = gamma11,
                       [4,7] = gamma12,
                       [4,8] = gamma13,
                       [4,9] = gamma14
        ;
        matrix _PHI_ [3,1] = 0,
                     [3,2] = 0,
                     [4,1] = 0,
                     [4,2] = 0,
                     [4,3] = 0,
                     [5,1] = 0,
                     [5,2] = 0,
                     [5,3] = 0,
                     [5,4] = 0,
                     [6,1] = 0,
                     [6,2] = 0,
                     [6,3] = 0,
                     [6,4] = 0,
                     [6,5] = 0,
                     [7,1] = 0,
                     [7,2] = 0,
                     [7,3] = 0,
                     [7,4] = 0,
                     [7,5] = 0,
                     [8,1] = 0,
                     [8,2] = 0,
                     [8,3] = 0,
                     [8,4] = 0,
                     [8,5] = 0,
                     [8,6] = 0,
                     [8,7] = 0,
                     [9,1] = 0,
                     [9,2] = 0,
                     [9,3] = 0,
                     [9,4] = 0,
                     [9,5] = 0,
                     [9,6] = 0,
                     [9,7] = 0,
                     [9,8] = 0
        ;
        run;
    proc printto PRINT=PRINT;
    run;
