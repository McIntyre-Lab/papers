/*******************************************************************************
* Filename: misc_tests_simultaneously_add_nonA.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I want to know what happens if I add all of the genes associated
* with nonA and heph
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

/* Splicing factor saturation with MMC */
    options ps = 70 ls=80;
    proc printto LOG=LOG PRINT='!MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_nona_nonc_heph.lst' NEW; run;
    proc calis data=SEM.cegsV_by_gene_sbs method=fiml COVARIANCE MEANSTR MODIFICATION maxiter=10000 outfit=fitstat;
        path
            fl_2_d Spf45 snf vir -> Sxl,
            Sxl -> FBgn0004227,
            FBgn0263968 FBgn0011224 -> Sxl,
            fl_2_d vir Sxl -> tra,
            tra2 tra -> fru dsx,
            dsx her ix -> Yp2 ,
            dsx -> FBgn0004367,

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
