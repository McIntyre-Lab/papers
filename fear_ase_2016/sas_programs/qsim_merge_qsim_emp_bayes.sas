/*******************************************************************************
* Filename: qsim_merge_qsim_emp_bayes.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I wan to compare qsim and emp bayes, so merge these results
* together
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname luis '!MCLAB/CEGS_Bayesian_Analysis_Paper/SAS_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge Qsim results with starting qsim value */
    proc sort data=CEGS.qsim_bayesian_results;
        by line fusion_id;
        run;

    proc sort data=CEGS.ase_qsim_tester;
        by line fusion_id;
        run;

    data qsim;
        merge CEGS.qsim_bayesian_results (in=in1) CEGS.ase_qsim_tester (in=in2);
        by line fusion_id;
        if in1;
        run;

/* Merge qsim and emp results */
    proc sort data=CEGS.emp_bayesian_results;
        by line mating_status fusion_id;
        run;

    proc sort data=qsim;
        by line mating_status fusion_id;
        run;

    proc sort data=CEGS.emp_bayesian_input;
        by line mating_status fusion_id;
        run;

    /* Create dataset and add flag_combined */
        * When qsim != 1 and qsim != 0 and abs(qsim - 0.5) > 0.05, then use qsim model
        * to estimate AI.
        *
        * Otherwise require all three empirical models combined to estimate AI.
        ;
        data CEGS.qsim_emp_theta_w_flag;
            retain line mating_status fusion_id q4_mean_theta q5_mean_theta
            q6_mean_theta qsim_mean_theta flag_q4_AI flag_q5_AI flag_q6_AI
            flag_all_AI flag_AI_qsim flag_AI_combined;
            merge CEGS.emp_bayesian_results (in=in1) qsim (in=in2) CEGS.emp_bayesian_input (in=in3);
            by line mating_status fusion_id;
            if in1;
            if qsim_tester ne 0 and qsim_tester ne 1 and abs(qsim_tester - 0.5) > 0.05 then do;
                if flag_AI_qsim eq 1 then flag_AI_combined = 1;
                else flag_AI_combined = 0;
            end;
            else do;
                if flag_q4_AI eq 1 and flag_q5_AI eq 1 and flag_q6_AI eq 1 then flag_AI_combined = 1;
                else flag_AI_combined = 0;
            end;
            drop q4_q: q5_q: q6_q: qsim_q: Bayesianpvalue_qsim LINE_TOTAL_: TESTER_TOTAL_: BOTH_TOTAL_: flag_analyze apn_:;
            run;

/* Export */
    proc export data=CEGS.qsim_emp_theta_w_flag outfile='!MCLAB/cegs_ase_paper/pipeline_output/ase_bayesian_qsim_emp_theta_w_flag.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
    proc datasets nolist;
        delete qsim;
        run; quit;
