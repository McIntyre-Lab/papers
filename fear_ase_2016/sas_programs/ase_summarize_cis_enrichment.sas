/*******************************************************************************
* Filename: ase_summarize_cis_enrichment.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: A set of enrichment tests to determine if there is an enrichment
* of cis-effects.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

%let flag=flag_AI_combined_m;

proc sort data=CEGS.clean_ase_sbs;
    by fusion_id;
    run;

proc freq data=CEGS.clean_ase_sbs noprint;
    by fusion_id;
    tables &flag / out=freq;
    run;

proc transpose data=freq out=flip;
    by fusion_id;
    var count;
    id &flag;
    run;

data clean;
    set flip;
    if _0 eq . then _0 = 0;
    if _1 eq . then _1 = 0;
    keep fusion_id _0 _1;
    run;

data flag_cis;
    set clean;
    if _1 > 0 then flag_cis_1 = 1;
    else flag_cis_1 = 0;
    if _1 >= 5 then flag_cis_5 = 1;
    else flag_cis_5 = 0;
    if _1 >= 10 then flag_cis_10 = 1;
    else flag_cis_10 = 0;
    if _1 >= 20 then flag_cis_20 = 1;
    else flag_cis_20 = 0;
    if _1 >= 30 then flag_cis_30 = 1;
    else flag_cis_30 = 0;
    if _1 >= 40 then flag_cis_40 = 1;
    else flag_cis_40 = 0;
    run;

ods rtf file='!MCLAB/cegs_ase_paper/pipeline_output/ase_summary/cis_enrichment.rtf';
proc freq data=flag_cis ;
    tables flag_cis_1 / chisq ;
    run;

proc freq data=flag_cis ;
    tables flag_cis_5 / chisq;
    run;

proc freq data=flag_cis ;
    tables flag_cis_10 / chisq;
    run;

proc freq data=flag_cis ;
    tables flag_cis_20 / chisq;
    run;

proc freq data=flag_cis ;
    tables flag_cis_30 / chisq;
    run;

proc freq data=flag_cis ;
    tables flag_cis_40 / chisq;
    run;
ods rtf close;

/* Merge to GO */
proc sort data=flag_cis;
    by fusion_id;
    run;

proc sort data=DMEL551.fusions2go;
    by fusion_id;
    run;

data merged;
    merge flag_cis (in=in1) DMEL551.fusions2go (in=in2);
    by fusion_id;
    if in1;
    run;

data CEGS.cis_line_flags_for_go;
    set merged;
    run;
