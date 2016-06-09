/********************************************************************************
* Import Bayesian summary and make flag
* 
* libname cegs '!MCLAB/cegs_sergey/sas_data';
* libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
* filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
* options SASAUTOS=(sasautos mymacros);
********************************************************************************/

 data WORK.BAYES    ;
 infile '!MCLAB/cegs_sergey/pipeline_output/ase_bayesian/bayesian_results_summary.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
    informat line $4. ;
    informat mating_status $1. ;
    informat fusion_id $9. ;
    informat q4_mean_theta best32. ;
    informat q4_q025 best32. ;
    informat q4_q975 best32. ;
    informat q5_mean_theta best32. ;
    informat q5_q025 best32. ;
    informat q5_q975 best32. ;
    informat q6_mean_theta best32. ;
    informat q6_q025 best32. ;
    informat q6_q975 best32. ;
    format line $4. ;
    format mating_status $1. ;
    format fusion_id $9. ;
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

data CEGS.flag_empirical_bayes;
    set bayes;
    if q4_q025 gt 0.5 or q4_q975 lt 0.5 then flag_q4_AI = 1;
    else flag_q4_AI = 0;
    if q5_q025 gt 0.5 or q5_q975 lt 0.5 then flag_q5_AI = 1;
    else flag_q5_AI = 0;
    if q6_q025 gt 0.5 or q6_q975 lt 0.5 then flag_q6_AI = 1;
    else flag_q6_AI = 0;
    run;
