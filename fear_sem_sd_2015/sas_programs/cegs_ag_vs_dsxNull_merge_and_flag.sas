/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge on gene2fbgn */
    proc sort data=SEM.cegsV_ag_yp2_added_genes;
        by gene;
        run;

    proc sort data=SEM.cegsV_gene2fbgn;
        by gene;
        run;

    data cegs_added;
        merge SEM.cegsV_ag_yp2_added_genes (in=in1) SEM.cegsV_gene2fbgn (in=in2);
        by gene;
        if in1;
        run;

/* merge dsxNull and cegs Adding genes */
    data dsx_repressed;
        set SEM.Dsxnullf_repressed_w_fb551;
        rename fb551 = primarY_fbgn;
        drop fb530;
        run;

    data dsx_induced;
        set SEM.Dsxnullf_induced_w_fb551;
        rename fb551 = primary_fbgn;
        drop fb530;
        run;

    proc sort data=dsx_repressed nodupkey;
        by primary_fbgn;
        run;

    proc sort data=dsx_induced nodupkey;
        by primary_fbgn;
        run;

    proc sort data=cegs_added;
        by primary_fbgn;
        run;

    data SEM.flag_ag_dsxNull;
        merge dsx_repressed (in=in1) dsx_induced (in=in2) cegs_added (in=in3);
        by primary_fbgn;
        if in3 then flag_add = 1;
        else flag_add = 0;
        if in1 then flag_dsxNull_repressed = 1;
        else flag_dsxNull_repressed =0;
        if in2 then flag_dsxNull_induced = 1;
        else flag_dsxNull_induced =0;
        run;

/* Export gene list of dsxNull repressed */
data tmp;
    set SEM.dsxnullf_repressed_w_fb551 ;
    drop fb530;
    run;

proc export data=tmp outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/dsxNullF_repressed_fb551.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean up */
proc datasets ;
    delete cegs_added;
    delete dsx_induced;
    delete dsx_repressed;
    delete tmp;
    run; quit;
