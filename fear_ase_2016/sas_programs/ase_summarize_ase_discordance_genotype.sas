/*******************************************************************************
* Filename: ase_summarize_ase_discordance_genotype.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Want to see if there are groups of genotypes that are biased
* towards the line and then others that are biased towards the tester.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data ai;
    set CEGS.clean_ase_sbs;
    if flag_AI_combined_m eq 1 and flag_AI_combined_v eq 1;
    if q5_mean_theta_m < 0.5 then ind_mated = 'mated_line';
    else if q5_mean_theta_m > 0.5 then ind_mated = 'mated_tester';
    if q5_mean_theta_v < 0.5 then ind_virgin = 'virgin_line';
    else if q5_mean_theta_v > 0.5 then ind_virgin = 'virgin_tester';
    drop flag_AI_combined_m flag_AI_combined_v;
    run;

proc sort data=ai; 
    by fusion_id;
    run;

proc freq data=ai noprint;
    by fusion_id;
    tables ind_mated /out=mfreq;
    tables ind_virgin /out=vfreq;
    run;

proc transpose data=mfreq out=mflip;
    by fusion_id;
    var count;
    id ind_mated;
    run;

data mated_geno_dis;
    set mflip;
    if mated_test eq . then delete;
    if mated_line eq . then delete;
    mated_total = mated_test + mated_line;
    drop _name_ _label_;
    run;

proc transpose data=vfreq out=vflip;
    by fusion_id;
    var count;
    id ind_virgin;
    run;

data virgin_geno_dis;
    set vflip;
    if virgin_test eq . then delete;
    if virgin_line eq . then delete;
    virgin_total = virgin_test + virgin_line;
    drop _name_ _label_;
    run;

data merged;
    merge mated_geno_dis (in=in1) virgin_geno_dis (in=in2);
    by fusion_id;
    run;

data CEGS.discordant_genotypes;
    set merged;
    run;

proc sort data=CEGS.discordant_genotypes;
    by DESCENDING virgin_total;
    run;


/* Clean up */
    proc datasets nolist;
        delete ai;
        delete mated_geno_dis;
        delete merged;
        delete mflip;
        delete mfreq;
        delete vflip;
        delete vfreq;
        delete virgin_geno_dis;
        run; quit;

