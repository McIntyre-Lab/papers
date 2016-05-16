/********************************************************************************
* Import Bayesian summary and make flag
* 
* libname cegs '!MCLAB/cegs_sergey/sas_data';
* libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
* filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
* options SASAUTOS=(sasautos mymacros);
********************************************************************************/

 data WORK.BAYES    ;
 infile '!MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/PG_model/PG_emp_bayesian_results.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
    informat line $5. ;
    informat mating_status $1. ;
    informat fusion_id $20. ;
    informat q4_mean_theta best32. ;
    informat q4_q025 best32. ;
    informat q4_q975 best32. ;
    informat q5_mean_theta best32. ;
    informat q5_q025 best32. ;
    informat q5_q975 best32. ;
    informat q6_mean_theta best32. ;
    informat q6_q025 best32. ;
    informat q6_q975 best32. ;
    format line $5. ;
    format mating_status $1. ;
    format fusion_id $20. ;
    format q4_mean_theta best12. ;
    format q4_q025 best12. ;
    format q4_q975 best12. ;
    format q5_mean_theta best12. ;
    format q5_q025 best12. ;
    format q5_q975 best12. ;
    format q6_mean_theta best12. ;
    format q6_q025 best12. ;
    format q6_q975 best12. ;
 input
             line $
             mating_status $
             fusion_id $
             q4_mean_theta 
             q4_q025 
             q4_q975 
             q5_mean_theta 
             q5_q025 
             q5_q975 
             q6_mean_theta 
             q6_q025 
             q6_q975 
 ;
 run;

data CEGS.emp_bayesian_results;
    set bayes;
    if q4_mean_theta eq . then flag_q4_AI = .;
    else if q4_q025 gt 0.5 or q4_q975 lt 0.5 then flag_q4_AI = 1;
    else flag_q4_AI = 0;
    if q5_mean_theta eq . then flag_q5_AI = .;
    else if q5_q025 gt 0.5 or q5_q975 lt 0.5 then flag_q5_AI = 1;
    else flag_q5_AI = 0;
    if q6_mean_theta eq . then flag_q6_AI = .;
    else if q6_q025 gt 0.5 or q6_q975 lt 0.5 then flag_q6_AI = 1;
    else flag_q6_AI = 0;
    if flag_q4_AI = 1 and flag_q5_AI = 1 and flag_q6_AI = 1 then flag_all_AI = 1;
    else flag_all_AI = 0;
    run;

proc datasets nolist;
    delete bayes;
    run;quit;
