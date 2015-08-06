/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

%let numsims=100;

*Combine all genes into a single dataset;
%macro combine_bic();
    %do i=1 %to &numsims; 
        libname addlink "!HOME/tmp/dspr_adding_links_simulation/&i./sas_data";
                
            data tmp;
                format model $14.;
                format bic best12.;
                format simulation best12.;
                set addlink.none;
                simulation = "&i.";
                drop gene;
                run;

            proc append data=tmp base=combined;
                run;
    %end;
%mend;

proc printto log="!HOME/tmp/dspr_adding_links_simulation/combine.log" new; run;
%combine_bic();
proc printto log=log new; run;

proc sort data=combined;
    by simulation bic ;
    run;

data best;
    set combined;
    by simulation ;
    if first.simulation then output;
    run;

proc sort data=best;
    by simulation;
    run;

proc freq data=best;
    by simulation;
    tables model;
    run;

/* Pull Baseline BIC for looking at diff */
data baseline;
    set combined;
    rename bic = baseline;
    if model eq 'Model_Baseline';
    drop model;
    run;

proc sort data=baseline;
    by simulation;
    run;

proc sort data=best;
    by simulation;
    run;

data merged_base;
    merge best (in=in1) baseline (in=in2);
    by simulation;
    diff = baseline - bic;
    run;

ods graphics on;
ods pdf file="!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/dspr_adding_links_simulation_bic_dist.pdf";
title "Distribution of difference of BIC from Baseline";
proc sgplot data=merged_base;
    density diff / type=kernel;
    run;

proc means data=merged_base min mean median max;
var diff;
run;
    
ods html;
proc freq data=merged_base;
table diff;
run;

ods pdf close;
ods graphics off;

/* Clean up */
proc datasets ;
    delete BASELINE;
    delete BEST;
    delete COMBINED;
    delete MERGED_BASE;
    delete TMP;
    run; quit;
