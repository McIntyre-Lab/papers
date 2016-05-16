/*******************************************************************************
* Filename: emp_bayesian_merge_gene_symbol_and_apn.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Look at the distribution of the thetas to determine how the
* empirical Bayesian estimates performed. Going to merge on APN and flags from
* the 100 genome simulation to remove fusions that will not perform well in the
* bayesian machine.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname MYCEGS '!HOME/storage/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge on apn information from input data */
    proc sort data=CEGS.emp_bayesian_results;
        by line mating_status fusion_id;
        run;

    proc sort data=CEGS.emp_bayesian_input;
        by line mating_status fusion_id;
        run;

    data merged;
        merge CEGS.emp_bayesian_results (in=in1) CEGS.emp_bayesian_input (in=in2);
        by line mating_status fusion_id;
        if in1;
        drop LINE_TOTAL_: TESTER_TOTAL_: apn_: flag_analyze;
        run;

    /* Create ranks based on apn */
        data rank_apn;
            set merged;
            if mean_apn < 5 then rank_apn = 0;
            else if mean_apn >= 5 and mean_apn < 15 then rank_apn = 1;
            else if mean_apn >= 15 and mean_apn < 50 then rank_apn = 2;
            else if mean_apn >= 50 then rank_apn = 3;
            run;

/* Merge on 100 Genome Simulation bias counts and make flag if a fusion showed bias */
    * NOTE: right now I have it set to if a fusion showed bias in any number of
    * simulated lines. I may want to tweak this to say 50% of lines;

    proc sort data=rank_apn;
        by fusion_id;
        run;

    proc sort data=CEGS.fb551_100_genome_bias_counts;
        by fusion_id;
        run;

    data merged2;
        merge rank_apn (in=in1) CEGS.fb551_100_genome_bias_counts (in=in2);
        by fusion_id;
        if in1;
        if lines_w_bias gt 0 then flag_fusion_biased = 1;
        else flag_fusion_biased = 0;
        drop lines_w_no_bias lines_w_no_allelic_reads lines_w_bias;
        run;

/* Merge on qsim and flag if qsim ne 0.5 */
    data qsim;
        set CEGS.ase_qsim_line;
        if qsim_line ne 0.5 then flag_qsim_bias = 1;
        else flag_qsim_bias = 0;
        run;

    proc sort data=qsim;
        by line fusion_id;
        run;

    proc sort data=merged2;
        by line fusion_id;
        run;

    data merged3;
        merge merged2 (in=in1) qsim (in=in2);
        by line fusion_id;
        if in1 then output merged3;
        run;

/* Create permanent dataset */
    data CEGS.emp_bayesian_results_w_flags;
        set merged3;
        run;

/* data for plot */
    data out;
        set CEGS.emp_bayesian_results_w_flags;
        rename q4_mean_theta = q4;
        rename q5_mean_theta = q5;
        rename q6_mean_theta = q6;
        rename mating_status = ms;
        *keep line mating_status fusion_id q4_mean_theta q5_mean_theta q6_mean_theta flag_q4_AI flag_q5_AI flag_q6_AI flag_all_AI rank_apn mean_apn flag_fusion_biased flag_qsim_bias qsim_line;
        run;


/* Export dataset */
    proc export data=out outfile='!MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/PG_model/emp_for_plotting.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Clean UP */
proc datasets nolist;
    delete merged;
    delete merged2;
    delete merged3;
    delete out;
    delete rank_apn;
    delete qsim;
    run; quit;

