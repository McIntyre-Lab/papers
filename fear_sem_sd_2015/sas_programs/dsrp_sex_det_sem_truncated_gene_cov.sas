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

/* SEM Gene Covariance */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/truncated_gene_cov_estimates.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_gene_cov_model.lst' NEW;
    run;
    proc calis data=w_sample method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d_A fl_2_d_B Spf45 snf vir -> Sxl_A,
            fl_2_d_A fl_2_d_B vir Sxl_A -> tra,
            tra2_A tra2_B tra -> fru,

            fl_2_d_A <-> snf  =0,
            fl_2_d_A <-> Spf45 =0,
            fl_2_d_A <-> tra2_A  =0,
            fl_2_d_A <-> tra2_B  =0,
            fl_2_d_A <-> vir =0,

            fl_2_d_B <-> snf  =0,
            fl_2_d_B <-> Spf45 =0,
            fl_2_d_B <-> tra2_A  =0,
            fl_2_d_B <-> tra2_B  =0,
            fl_2_d_B <-> vir =0,

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

    proc export data=combined outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_gene_cov_estimates.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Create perminant dataset for BIC */
    * I need this later so this was the best place to output it;
    data SEM.dsrp_sex_det_tgene_cov_model_bic;
        length model $12.;
        retain model;
        set fitstat;
        where IndexCode = 312;
        model = "Baseline";
        rename FitValue = BIC;
        label FitValue = ' ';
        keep model FitValue;
        run;

/* Clean up */
    proc datasets nolist;
        delete combined;
        delete covstd;
        delete pathstd;
        delete pathstd2;
        delete w_sample;
        delete fitstat;
        run; quit;
