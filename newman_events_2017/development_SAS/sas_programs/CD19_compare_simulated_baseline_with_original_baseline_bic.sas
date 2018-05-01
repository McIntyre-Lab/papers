/* Libraries */

ods html close; ods listing;
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname con '!PATCON/sas_data';
libname sem '!MCLAB/grants/ase_grn/INSR-mTOR_T1D/sas_data';

/* Combine results of simulations to compare BICs with original path */

%let numsims=100;

*Combine all genes into a single dataset;
%macro combine_bic();
    %do i=1 %to &numsims; 
        libname addgen "/mnt/store/grn_sandbox/cd19_8000_baseline_simulation/8000_gene_sim/&i./sas_data";

            data tmp;
                format gene $36.;
                format model $14.;
                format bic best12.;
                format simulation best12.;
                format numgenes best12.;
                set addgen.cd19;
                simulation = "&i.";
                numgenes = "&numgenes.";
                run;

            proc append data=tmp base=combined;
                run;
    %end;
%mend;

proc printto log="/mnt/store/grn_sandbox/CD19_100_sim_combine.log" new; run;
%combine_bic();
proc printto log=log new; run;

proc sort data=combined;
    by simulation gene bic ;
    run;

/* Add in original baseline BIC and subtract from simulated baseline */


data diff;
   set combined;
   genecov_bic=1711.4367;
   fullcov_bic=1800.8199;
   genecov_diff=genecov_bic - bic;
   fullcov_diff=fullcov_bic - bic;
run;

ods graphics on;
ods pdf file="!MCLAB/grants/ase_grn/INSR-mTOR_T1D/analysis_output/CD19_100simulations_Baseline_only_sim_bic_dist_fullcovar.pdf";
title "Distribution of difference of simulated Baseline BIC from gene covariance model Baseline";
proc sgplot data=diff;
    density genecov_diff / type=kernel;
    run;

proc means data=diff min mean median max;
var genecov_diff;
run;
    
ods html;
proc freq data=diff;
table genecov_diff;
run;

ods pdf close;
ods graphics off;



ods graphics on;
ods pdf file="!MCLAB/grants/ase_grn/INSR-mTOR_T1D/analysis_output/CD19_100simulations_Baseline_only_sim_bic_dist_fullcovar.pdf";
title "Distribution of difference of simulated Baseline BIC from full covariance model Baseline";
proc sgplot data=diff;
    density fullcov_diff / type=kernel;
    run;

proc means data=diff min mean median max;
var fullcov_diff;
run;
    
ods html;
proc freq data=diff;
table fullcov_diff;
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

