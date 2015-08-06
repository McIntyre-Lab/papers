/*
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel530 "!MCLAB/useful_dmel_data/flybase530/sas_data";
*/


/* Sort Datasets by FBgn_cat */

data FBGn;
    set dmel530.Fb530_si_fusions_unique_flagged;
    keep fbgn_cat;
    run;

proc sort data=FBGN nodupkey;
    by fbgn_cat;
    run;

proc sort data=SEM.dsxD_induced;
    by fbgn_cat;
    run;

proc sort data=SEM.dsxD_repressed;
    by fbgn_cat;
    run;

proc sort data=SEM.dsxNullF_induced;
    by fbgn_cat;
    run;

proc sort data=SEM.dsxNullF_repressed;
    by fbgn_cat;
    run;

proc sort data=SEM.dsx_ctrl_female_induced;
    by fbgn_cat;
    run;

proc sort data=SEM.dsx_ctrl_female_repressed;
    by fbgn_cat;
    run;

data merged;
    merge FBGN (in=in1) SEM.dsxD_induced (in=in2) SEM.dsxD_repressed (in=in3) 
    SEM.dsxNullF_induced (in=in4) SEM.dsxNullF_repressed (in=in5) 
    SEM.dsx_ctrl_female_induced (in=in6) SEM.dsx_ctrl_female_repressed (in=in7);
    by fbgn_cat;
    if in1 then do;
        if in2 then flag_dsxd = 1; 
        else if in3 then flag_dsxd = -1; 
        else flag_dsxd = 0;

        if in4 then flag_dsxNull = 1; 
        else if in5 then flag_dsxNull = -1; 
        else flag_dsxNull = 0;

        if in6 then flag_ctrl = 1; 
        else if in7 then flag_ctrl = -1; 
        else flag_ctrl = 0;
    end;
    run;

proc export data=merged outfile='/tmp/dsx_flags.csv' dbms=csv replace;
    putnames=yes;
    run;

data _null_;
    call system('Rscript $MCLAB/cegs_sem_sd_paper/r_programs/dsx_venn.R /tmp/dsx_flags.csv');
    run;
