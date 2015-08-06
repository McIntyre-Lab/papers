/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* CEGS ag DSX list */
    data ag;
        set SEM.cegsV_ag_flag_dsx;
        if flag_all_dsx_m4 = 1 or flag_ll_dsx_m26 = 1;
        rename primary_fbgn = fbgn;
        keep symbol primary_fbgn;
        run;

/* merge dsxNull and CEGS adding genes downstream of DSX */
    data dsx;
        set SEM.Dsxnullf_repressed_w_fb551;
        rename fb551 = fbgn;
        run;

    proc sort data=ag;
        by fbgn;
        run;

    proc sort data=dsx;
        by fbgn;
        run;

    data SEM.flag_dsxNull_ag_ds_dsx;
        merge dsx (in=in1) ag (in=in2);
        by fbgn;
        if in1 then flag_dsxnull = 1;
        else flag_dsxnull = 0;
        if in2 then flag_ds_dsx = 1;
        else flag_ds_dsx = 0;
        drop fb530
        run;

    proc freq data=SEM.flag_dsxNull_ag_ds_dsx;
        tables flag_dsxNull*flag_ds_dsx /chisq;
        run;

    data valid;
        set SEM.flag_dsxNull_ag_ds_dsx;
        if flag_dsxnull = 1 and flag_ds_dsx = 1;
        run;

