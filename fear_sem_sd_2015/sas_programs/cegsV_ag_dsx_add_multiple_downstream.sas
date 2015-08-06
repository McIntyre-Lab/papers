/*******************************************************************************
* Filename: cegsV_ag_dsx_add_multiple_downstream.sas
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
        set SEM.cegsV_ag_yp2_flag_ds_dsx;
        where flag_best_dsx_m3 = 1;
        run; * 13 genes;

    /* 13 gene are:
        symbol              primary_fbgn
        --------------------------------
        lab                  FBgn0002522
        sbr                  FBgn0003321
        mei-41               FBgn0004367
        RpIII128             FBgn0004463
        Hmgs                 FBgn0010611
        CG7099               FBgn0032517
        CG15909              FBgn0033090
        CG4080               FBgn0035983
        ETHR                 FBgn0038874
        CG10669              FBgn0039329
        v(2)k05816           FBgn0042627
        Snap                 FBgn0250791
        mxc                  FBgn0260789
    */

/* Baseline Check */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
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
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2
        ;
        run; * BIC = 455.8954;

/* Add All DS DSX */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0002522 FBgn0003321 FBgn0004367 FBgn0004463 FBgn0010611
            FBgn0032517 FBgn0033090 FBgn0035983 FBgn0038874 FBgn0039329 FBgn0042627
            FBgn0250791 FBgn0260789,

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
        run; * BIC = 303.9164;

/* Add All DS DSX fullcov model */
    proc calis data=SEM.cegsv_by_gene_sbs method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0002522 FBgn0003321 FBgn0004367 FBgn0004463 FBgn0010611
            FBgn0032517 FBgn0033090 FBgn0035983 FBgn0038874 FBgn0039329 FBgn0042627
            FBgn0250791 FBgn0260789
        ;
        run; * BIC = 178.9;

/* Compare 13 to DSX null list*/
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
        run;

/* Create gene list of genes added DS of DSX */
    data m325;
        set SEM.cegsV_ag_yp2_flag_ds_dsx;
        where flag_all_dsx_m25 = 1;
        run; * 31 genes;

    proc sort data=m325;
        by primary_fbgn;
        run;

    data merged;
        merge dsxNull(in=in1) m325 (in=in2);
        by primary_fbgn;
        if in1 and in2;
        run;
