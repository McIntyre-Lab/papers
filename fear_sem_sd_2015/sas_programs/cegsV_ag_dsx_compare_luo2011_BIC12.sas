/*******************************************************************************
* Filename: cegsV_ag_dsx_compare_luo2011_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Luo 2011 identified 58 genes with dsx binding sites in them. I
* want to compare their list to the genes added downstream of dsx.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge on DSX flags */
    proc sort data=SEM.luo2011;
        by primary_fbgn;
        run;

    proc sort data=SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12;
        by primary_fbgn;
        run;

    data merged2;
        merge SEM.luo2011 (in=in1) SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12 (in=in2);
        by primary_fbgn;
        if in1 then flag_luo = 1;
        else flag_luo = 0;
        run;

/* Do the Freqs */
    proc freq data=merged2;
        tables flag_luo*flag_all_dsx_m3;
        tables flag_luo*flag_all_dsx_m25;
        tables flag_luo*flag_best_dsx_m3;
        tables flag_luo*flag_best_dsx_m25;
        run;

    data peek;
        set merged2;
        if flag_luo = 1 and (flag_all_dsx_m3 = 1 or flag_all_dsx_m25 =1);
        run;

    proc print data=peek;
    run;

/* 
                   primary_                   all_dsx_    all_dsx_    dsx_model3_    dsx_model25_     best_     best_
Obs    symbol         fbgn           flag_sig       m3          m25          rank           rank        dsx_m3    dsx_m25    flag_luo
1     ds             FBgn0000497        1           1           0             9              .            0         0           1
2     GEFmeso        FBgn0050115        0           1           1            28             34            0         0           1

*/


data peek2;;
    set merged2;
    if flag_all_dsx_m3 = 1 or flag_all_dsx_m25 = 1 then flag_dsx = 1;
    else flag_dsx = 0;
    run;

proc freq data=peek2;
    tables flag_dsx*flag_luo;
    run;


/* Clean up */
proc datasets ;
    delete luo;
    delete merged;
    delete merged2;
    delete oops;
    delete peek;
    delete peek2;
    run; quit;


