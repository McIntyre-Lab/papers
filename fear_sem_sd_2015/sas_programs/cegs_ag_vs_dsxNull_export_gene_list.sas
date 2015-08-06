/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* How many genes overlap between the different sets */
    proc freq data=SEM.flag_ag_dsxNull;
        tables flag_dsxNull_repressed*flag_add; * 26 overlap;
        tables flag_dsxNull_induced*flag_add; * 31 overlap;
        run;

    data repressed;
        set SEM.flag_ag_dsxNull;
        if flag_add =1 and flag_dsxNull_repressed = 1;
        keep primary_fbgn model;
        run;

    data induced;
        set SEM.flag_ag_dsxNull;
        if flag_add =1 and flag_dsxNull_induced = 1;
        keep primary_fbgn model;
        run;

/* merge repressed and gene symbol */
    proc sort data=repressed;
        by primary_fbgn;
        run;

    proc sort data=DMEL551.symbol2fbgn;
        by primary_fbgn;
        run;

    data merged;
        retain symbol primary_fbgn;
        merge repressed (in=in1) DMEL551.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1;
        run;

    proc export data=merged outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_overlap_w_dsxNull_repressed.csv' dbms=csv replace;
        putnames=yes;
        run;

/* merge induced and gene symbol */
    proc sort data=induced;
        by primary_fbgn;
        run;

    proc sort data=DMEL551.symbol2fbgn;
        by primary_fbgn;
        run;

    data merged;
        retain symbol primary_fbgn;
        merge induced (in=in1) DMEL551.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1;
        run;

    proc export data=merged outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes_vs_dsxNull/cegs_ag_overlap_w_dsxNull_induced.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
    proc datasets;
        delete merged;
        delete induced;
        delete repressed;
        run; quit;
