/*******************************************************************************
* Filename: cegsV_ag_dsx_compare_luo2011.sas
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

    proc sort data=SEM.cegsV_ag_yp2_flag_ds_dsx;
        by primary_fbgn;
        run;

    data merged2;
        merge SEM.luo2011 (in=in1) SEM.cegsV_ag_yp2_flag_ds_dsx (in=in2);
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

    data peak;
        set merged2;
        if flag_luo = 1 and (flag_all_dsx_m3 = 1 or flag_all_dsx_m25 =1);
        run;

    proc print data=peak;
    run;

/* Clean up */
proc datasets ;
    delete luo;
    delete merged;
    delete merged2;
    delete oops;
    delete peak;
    run; quit;
