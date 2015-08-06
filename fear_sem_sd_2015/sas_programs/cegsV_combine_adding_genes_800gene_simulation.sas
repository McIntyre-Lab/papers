/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

%let numgenes=800;
%let numsims=10;

*Combine all genes into a single dataset;
%macro combine_bic();
    %do i=1 %to &numsims; 
        libname addgen "!HOME/tmp/cegs_adding_genes_simulation/&numgenes._gene_sim/&i./sas_data";
        %do j=1 %to &numgenes;
                
            data tmp;
                format gene $15.;
                format model $14.;
                format bic best12.;
                format simulation best12.;
                format numgenes best12.;
                set addgen.gene&j;
                simulation = "&i.";
                numgenes = "&numgenes.";
                run;
            proc append data=tmp base=combined;
                run;
        %end;
    %end;
%mend;

proc printto log="!HOME/tmp/cegs_adding_genes_simulation/&numgenes._gene_sim/combine.log" new; run;
%combine_bic();
proc printto log=log new; run;

proc sort data=combined;
    by simulation gene bic ;
    run;

data best;
    set combined;
    by simulation gene;
    if first.gene then output;
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
    drop model numgenes;
    run;

proc sort data=baseline;
    by gene simulation;
    run;

proc sort data=best;
    by gene simulation;
    run;

data merged_base;
    merge best (in=in1) baseline (in=in2);
    by gene simulation;
    diff = baseline - bic;
    run;

ods graphics on;
ods pdf file="!MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegs_adding_genes_simulation_&numgenes._gene_sim_bic_dist.pdf";
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
    delete REGSTRY;
    delete SASMACR;
    delete TMP;
    run; quit;

