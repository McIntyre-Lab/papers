/*******************************************************************************
* Filename: tier_misspecification_table.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: This script calculates TIER and formats a table for comparing
* emp Bayes vs Binomial.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab needed columns */
    data miss;
        retain FUSION_ID q_true;
        set CEGS.r101_misspecification;
        if flag_AI_q4 = 1 and flag_AI_q5 = 1 and flag_AI_q6 = 1 then flag_AI_emp = 1;
        else flag_AI_emp = 0;
        keep FUSION_ID q_true flag_AI_binomial_test flag_AI_q4 flag_AI_q5 flag_AI_q6 flag_AI_emp;
        run;

/* Calculate the TIER for each one */
    proc freq data=miss noprint;
        by q_true;
        table flag_AI_binomial_test /out=bin;
        table flag_AI_q4 /out=q4;
        table flag_AI_q5 /out=q5;
        table flag_AI_q6 /out=q6;
        table flag_AI_emp /out=qemp;
        run;

/* Clean up FREQ output to just have TIER */
    %macro clean(mydat, flag, name);
        data &mydat._2;
            set &mydat;
            where &flag eq 1 ;
            &name = percent / 100;
            per_miss = ((q_true / 0.5) - 1) * 100;
            drop &flag count percent q_true;
            run;

        proc datasets nolist;
            delete &mydat;
            run; quit;

    %mend;
    %clean(bin,flag_AI_binomial_test, binomial);
    %clean(q4,flag_AI_q4, pg4);
    %clean(q5,flag_AI_q5, pg5);
    %clean(q6,flag_AI_q6, pg6);
    %clean(qemp,flag_AI_emp, pgemp);

/* Merge all of the different levels of bias together */
    data miss_bias;
        retain per_miss;
        merge bin_2 (in=in1) q4_2 (in=in2) q5_2 (in=in3) q6_2 (in=in4) qemp_2 (in=in5);
        by per_miss;
        run;

/* Export Dataset for plotting */
proc export data=miss_bias outfile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/r101_misspecification_tier_binomial_vs_pg_plotting.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean Up */
proc datasets nolist;
    delete bin_2;
    delete miss;
    delete miss_bias;
    delete q4_2;
    delete q5_2;
    delete q6_2;
    delete qemp_2;
    run; quit;


