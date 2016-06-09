/*******************************************************************************
* Filename: qsim_filter_luis_qsim_qdna_results.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: In the Luis paper we filtered the results to keep fusions where 
*
* qSIM != 0
* qSIM != 1 
* abs(qSIM - .5) > 0.2
* 
* In the paper we say that n = 617, but I find here n = 317. I am
* thinking that 617 was a typo.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname luis '!MCLAB/CEGS_Bayesian_Analysis_Paper/SAS_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

data fusions_w_qsim;
    set luis.berlin_c1674_new_algn_simulation;
    where p_new ne 0 and p_new ne 1 and p_new ne .;
    if abs(p_new - 0.5) > 0.1;
    keep fusion_id;
    run;

proc sort data=fusions_w_qsim;
    by fusion_id;
    run;

proc sort data=CEGS.luis_pg_results;
    by fusion_id;
    run;

data merged;
    merge CEGS.luis_pg_results (in=in1) fusions_w_qsim (in=in2);
    by fusion_id;
    if in1 and in2;
    if mean_PG_q_DNA ne . and mean_PG_q_sim ne . and mean__PG_q_ahalf ne .;
    rename mean__PG_q_ahalf = mean_PG_q_ahalf;
    keep fusion_id mean_PG_q_DNA mean_PG_q_sim mean__PG_q_ahalf;
    run;

proc export data=merged outfile='!MCLAB/cegs_ase_paper/pipeline_output/qsim/output/luis_paper_pg_qsim_vs_qdna.csv' dbms=csv replace;
    putnames=yes;
    run;

