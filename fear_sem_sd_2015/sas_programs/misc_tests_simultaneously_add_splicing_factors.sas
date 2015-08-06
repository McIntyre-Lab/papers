/*******************************************************************************
* Filename: misc_tests_simultaneously_add_splicing_factors.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I want to know if I can reach saturation to the F1-hybrid
* dataset. One of our arguments why DSPR did not add any genes is that there
* was not enough variation, therefore the model was saturated. I am going to
* try adding all of the splicing factors simultaneously to the model.
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Baseline model */
    options ps = 70 ls=80;
    proc printto LOG='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/baseline_model.log' PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/baseline_model.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
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
        run;

/* Identify set of splicing factors to add */
    options ps = 100 ls=200;
    proc printto print='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_corr_analysis.lst' new; run;
    proc corr data=SEM.cegsV_by_gene_sbs out=corr;
        var FBgn0003261 FBgn0031883 FBgn0019990 FBgn0039600 FBgn0037660 FBgn0037466
        FBgn0027548 FBgn0014366 FBgn0050122 FBgn0011224 FBgn0028474 FBgn0004227
        FBgn0261641 FBgn0038344 FBgn0263603 FBgn0005649 FBgn0030625 FBgn0000146
        FBgn0027587 FBgn0022942 FBgn0036734 FBgn0036915 FBgn0035253 FBgn0016978
        FBgn0036386;
        run;

    proc printto print='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_factor_analysis.lst' new; run;
    proc factor data=SEM.cegsV_by_gene_sbs simple method=PRINCIPAL scree rotate=varimax nfact=10 round flag=.4 outstat=factor;
                var FBgn0003261 FBgn0031883 FBgn0019990 FBgn0039600 FBgn0037660
                FBgn0037466 FBgn0027548 FBgn0014366 FBgn0050122 FBgn0011224
                FBgn0028474 FBgn0004227 FBgn0261641 FBgn0038344 FBgn0263603
                FBgn0005649 FBgn0030625 FBgn0000146 FBgn0027587 FBgn0022942
                FBgn0036734 FBgn0036915 FBgn0035253 FBgn0016978 FBgn0036386;
                run;

    * Export for MMC;
    proc transpose data=SEM.cegsV_by_gene_sbs out=flip;
        var FBgn0003261 FBgn0031883 FBgn0019990 FBgn0039600 FBgn0037660
        FBgn0037466 FBgn0027548 FBgn0014366 FBgn0050122 FBgn0011224 FBgn0028474
        FBgn0004227 FBgn0261641 FBgn0038344 FBgn0263603 FBgn0005649 FBgn0030625
        FBgn0000146 FBgn0027587 FBgn0022942 FBgn0036734 FBgn0036915 FBgn0035253
        FBgn0016978 FBgn0036386;
        id line;
        run;

    proc export data=flip outfile='/home/jfear/sandbox/splice.csv' dbms=csv replace;
        putnames=yes;
        run;

    /* MMC module groupings

    FBgn0011224 FBgn0004227
            
    FBgn0038344 FBgn0030625

    FBgn0035253 FBgn0261641 FBgn0039600 FBgn0036734 FBgn0014366 FBgn0028474 FBgn0027587 FBgn0016978 FBgn0037660

    FBgn0050122 FBgn0019990 FBgn0003261 FBgn0000146 FBgn0005649 FBgn0022942 FBgn0036386

    FBgn0263603 FBgn0037466 FBgn0031883 FBgn0027548 FBgn0036915

    */

/* Splicing factor saturation with MMC */
    options ps = 70 ls=80;
    proc printto LOG=LOG PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_mmc.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            Sxl -> FBgn0011224  FBgn0038344 FBgn0035253 FBgn0050122 FBgn0263603,
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
        run;

/* Splicing factor saturation with Factor 1 */
    options ps = 70 ls=80;
    proc printto LOG=LOG PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_factor.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            Sxl -> FBgn0003261 FBgn0019990 FBgn0039600 FBgn0037466 FBgn0027548 FBgn0014366 FBgn0050122 FBgn0011224 FBgn0028474 FBgn0004227 FBgn0261641 FBgn0263603 FBgn0005649 FBgn0000146 FBgn0027587 FBgn0036734 FBgn0035253,
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
        run;

/* Identify set of genes that were added downstream of dsx */
    proc printto PRINT=PRINT new; run;

    proc print data=SEM.cegsV_ag_yp2_flag_ds_dsx_bic12 (where=(flag_best_dsx_m3 = 1 or flag_best_dsx_m25 = 1));
        run;

    proc transpose data=SEM.cegsV_by_gene_sbs out=dsxflip;
        var FBgn0002522 FBgn0004367 FBgn0019925 FBgn0028480 FBgn0032517 FBgn0032886 FBgn0034432 FBgn0035771 FBgn0039767 FBgn0040281 FBgn0250791 FBgn0260789;
        id line;
        run;

    * Export for MMC;
    proc export data=dsxflip outfile='/home/jfear/sandbox/dsx.csv' dbms=csv replace;
        putnames=yes;
        run;

    /* MMC module groupings

    FBgn0004367 FBgn0260789 FBgn0032517

    FBgn0034432 FBgn0028480 FBgn0019925 FBgn0035771 FBgn0039767

    FBgn0040281 FBgn0250791 FBgn0032886 FBgn0002522

    */

/* DSX saturation with MMC */
    proc printto LOG=LOG PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_dsx_mmc.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0004367 FBgn0034432 FBgn0040281 ,

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
        run;

/* Splicing factors and DSX saturation with MMC */
    options ps = 70 ls=80;
    proc printto LOG=LOG PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_or_dsx_mmc.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            Sxl -> FBgn0011224  FBgn0038344 FBgn0035253 FBgn0050122 FBgn0263603 FBgn0052473 FBgn0052702 ,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0004367 FBgn0034432 FBgn0040281 FBgn0035771 FBgn0002522 ,

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
        run;
