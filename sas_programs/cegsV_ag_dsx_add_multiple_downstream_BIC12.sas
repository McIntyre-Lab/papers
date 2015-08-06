/*******************************************************************************
* Filename: cegsV_ag_dsx_add_multiple_downstream_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I am wanting to see what happens when I try to add multiple
* genes simultaneously downstream of DSX, so I need to create a list of genes
* that had their best fit downstream of DSX (model 3: Dsx -> gene).
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create gene list of genes added DS of DSX */
    data m3;
        set SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12;
        where flag_best_dsx_m3 = 1;
        run; * 5 genes;

    /* 5 gene are:

        proc print data=m3;
        run;

        symbol              primary_fbgn
        --------------------------------
        lab                 FBgn0002522
        mei-41              FBgn0004367
        CG7099              FBgn0032517
        Snap                FBgn0250791
        mxc                 FBgn0260789
    */

/* Baseline Check */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE MEANSTR RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2  =0,
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
            snf <-> tra2  =0,
            snf <-> vir =0,

            Spf45 <-> tra2  =0,
            Spf45 <-> vir =0,

            tra2 <-> vir =0
        ;
        run; * BIC = 524.7484;

/* Full cov check */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE MEANSTR RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2
        ;
        run; * BIC = 507.6693;


/* Add All DS DSX */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE MEANSTR RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0002522 FBgn0004367 FBgn0032517 FBgn0250791 FBgn0260789,

            fl_2_d <-> her  =0,
            fl_2_d <-> ix  =0,
            fl_2_d <-> snf  =0,
            fl_2_d <-> Spf45 =0,
            fl_2_d <-> tra2  =0,
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
            snf <-> tra2  =0,
            snf <-> vir =0,

            Spf45 <-> tra2  =0,
            Spf45 <-> vir =0,

            tra2 <-> vir =0
        ;
        run; * BIC = 364.2416;

/* Add All DS DSX fullcov model */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE MEANSTR RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0002522 FBgn0004367 FBgn0032517 FBgn0250791 FBgn0260789
        ;
        run; * BIC = 347.1625;

/* Compare 5 to DSX null list*/
    data dsxNull;
        set SEM.Dsxnullf_repressed SEM.Dsxnullf_induced;
        rename fbgn_cat = primary_fbgn;
        run;

    proc sort data=dsxNull;
        by primary_fbgn;
        run;

    proc sort data=m3;
        by primary_fbgn;
        run;

    data merged;
        merge dsxNull(in=in1) m3 (in=in2);
        by primary_fbgn;
        if in1 and in2;
        run; * only CG7099 was added to model 3 and in null list;

/* Create gene list of genes added DS of DSX */
    data m325;
        set SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12;
        where flag_all_dsx_m3 = 1 or flag_all_dsx_m25 = 1;
        run; * 8 genes;

    proc sort data=m325;
        by primary_fbgn;
        run;

    data merged;
        merge dsxNull(in=in1) m325 (in=in2);
        by primary_fbgn;
        if in1 and in2;
        run; * 12 genes added downstream of dsx and in dsx null;

    /* 
    proc print data=merged;
    run;
                                 flag_       flag_                                     flag_     flag_
    primary_                      all_dsx_    all_dsx_    dsx_model3_    dsx_model25_     best_     best_
    fbgn           symbol            m3          m25          rank           rank        dsx_m3    dsx_m25

    FBgn0003204    ras                1           1             2              9            0         0
    FBgn0014184    Oda                1           1            31             34            0         0
    FBgn0020503    CLIP-190           1           1            22             36            0         0
    FBgn0024366    CG11409            1           1             4              8            0         0
    FBgn0030608    Lsd-2              0           1             .              2            0         0
    FBgn0032517    CG7099             1           1             1              8            1         0
    FBgn0035850    Atg18              1           1            17             33            0         0
    FBgn0036299    Tsf2               1           1            31             34            0         0
    FBgn0038659    endoA              0           1             .             17            0         0
    FBgn0050115    GEFmeso            1           1            28             34            0         0
    FBgn0260936    scny               1           1            35             36            0         0
    FBgn0261564    CG42678            1           1             4              7            0         0

    */
