/********************************************************************************
* libname cegs '!MCLAB/cegs_sergey/sas_data';
* libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
* libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
* filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
* options SASAUTOS=(sasautos mymacros);
********************************************************************************/

/* Merge on Gene Symbol */
proc sort data=dmel.Fb551_si_fusions_unique_flagged;
    by fusion_id;
    run;

proc sort data=CEGS.flag_empirical_bayes;
    by fusion_id;
    run;
    
data CEGS.flag_empirical_bayes_w_symbol;
    merge CEGS.flag_empirical_bayes (in=in1) dmel.Fb551_si_fusions_unique_flagged (in=in2);
    by fusion_id;
    if in1;
    keep line mating_status fusion_id q4_mean_theta q4_q025 q4_q975
    q5_mean_theta q5_q025 q5_q975 q6_mean_theta q6_q025 q6_q975 flag_q4_AI
    flag_q5_AI flag_q6_AI Genes_per_fusion symbol_cat FBgn_cat Constitutive
    Common Alternative ;
    run;

/* Export dataset */
proc export data=CEGS.flag_empirical_bayes_w_symbol outfile='!MCLAB/cegs_sergey/reports/ase/ase_bayesian_preliminary_results.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Create distribution plots of theta */
data _null_;
    call system('Rscript $MCLAB/cegs_sergey/r_programs/ase_bayesian_plot_mini_empirical_thetas.R');
    run;
