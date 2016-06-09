/*******************************************************************************
* Filename: tier_r101_misspecification_bias.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: This script calculates the TIER and formats the table for
* plotting.
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
        set CEGS.r101_fb557_miss;
        if flag_AI_q4 = 1 and flag_AI_q5 = 1 and flag_AI_q6 = 1 then flag_AI_emp = 1;
        else flag_AI_emp = 0;
        keep FUSION_ID q_true flag_AI_q_true flag_AI_q_sim_above_01per
        flag_AI_q_sim_above_02per flag_AI_q_sim_above_05per
        flag_AI_q_sim_above_10per flag_AI_q_sim_above_20per flag_AI_emp;
        run;

/* Calculate the TIER for each one */
    proc freq data=miss noprint;
        by q_true;
        table flag_AI_q_true /out=qtrue;
        table flag_AI_q_sim_above_01per /out=above01;
        table flag_AI_q_sim_above_02per /out=above02;
        table flag_AI_q_sim_above_05per /out=above05;
        table flag_AI_q_sim_above_10per /out=above10;
        table flag_AI_q_sim_above_20per /out=above20;
        table flag_AI_emp /out=qemp;
        run;

/* Clean up FREQ output to just have TIER */
    %macro clean(mydat, flag, name);
        data &mydat._2;
            set &mydat;
            where &flag eq 1;
            &name = percent / 100;
            sim_bias = 200 * (q_true - 0.5);
            drop &flag count percent q_true;
            run;

        proc datasets nolist;
            delete &mydat;
            run; quit;

    %mend;
    %clean(qtrue,flag_AI_q_true, tier_00);
    %clean(above01,flag_AI_q_sim_above_01per,tier_01);
    %clean(above02,flag_AI_q_sim_above_02per,tier_02);
    %clean(above05,flag_AI_q_sim_above_05per,tier_05);
    %clean(above10,flag_AI_q_sim_above_10per,tier_10);
    %clean(above20,flag_AI_q_sim_above_20per,tier_20);
    %clean(qemp,flag_AI_emp, tier_emp);

/* Merge all of the different levels of bias together */
    data miss_bias;
        retain sim_bias;
        merge 
        qtrue_2 (in=in1) 
        above01_2 (in=in2) 
        above02_2 (in=in3) 
        above05_2 (in=in4) 
        above10_2 (in=in5) 
        above20_2 (in=in6) 
        qemp_2 (in=in7) ;
        by sim_bias;
        run;

/* Export Dataset for plotting */
proc export data=miss_bias outfile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/r101_fb557_miss_bias_for_plotting.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean Up */
proc datasets nolist;
    delete above01_2;
    delete above02_2;
    delete above05_2;
    delete above10_2;
    delete above20_2;
    delete miss;
    delete miss_bias;
    delete qtrue_2;
    delete qemp_2;
    run; quit;
