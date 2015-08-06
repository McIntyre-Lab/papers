/*******************************************************************************
* Filename: cegsV_ag_yp2_vs_dsxNull_compare_luo2011.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Luo 2011 identified 58 genes with dsx binding sites in them. I
* want to compare their list to the genes added to the GRN.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge on Luo flags */
    proc sort data=SEM.luo2011;
        by primary_fbgn;
        run;

    proc sort data=SEM.flag_ag_dsxNull;
        by primary_fbgn;
        run;

    proc sort data=SEM.cegsV_ag_yp2_flag_ds_dsx;
        by primary_fbgn;

    data merged;
        retain primary_fbgn symbol model flag_dsxNull_repressed flag_dsxNull_induced flag_luo flag_sig;
        merge SEM.luo2011 (in=in1) SEM.flag_ag_dsxNull (in=in2) SEM.cegsV_ag_yp2_flag_ds_dsx (in=in3);
        by primary_fbgn;
        if in1 then flag_luo = 1;
        else flag_luo = 0;
        if flag_add = 1;
        rename flag_sig = flag_luo_damid_sig;
        drop gene;
        run; 

    proc export data=merged outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_w_flag_dsxNull_flag_luo.csv' dbms=csv replace;
        putnames=yes;
        run;

    data SEM.cegsv_ag_w_flag_ds_dsx_flag_luo;
        set merged;
        run;

data dsx;
    set SEM.cegsv_ag_w_flag_ds_dsx_flag_luo;
    if flag_dsxNull_repressed = 1 then output;
    if (dsx_model3_rank ne '.' and dsx_model3_rank <5) or (dsx_model25_rank ne '.' and dsx_model25_rank <5);
    run;

proc sort data=dsx;
    by dsx_model3_rank;
    run;



/* Clean up */
proc datasets ;
    delete merged;
    delete dsx;
    run; quit;
