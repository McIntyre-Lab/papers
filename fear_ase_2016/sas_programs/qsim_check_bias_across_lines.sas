/*******************************************************************************
* Filename: qsim_check_bias_across_lines.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Do the same genomic regions show bias across all lines.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Make wide dataset */
    proc sort data=CEGS.ase_qsim_tester;
        by fusion_id;
        run;

    proc transpose data=CEGS.ase_qsim_tester out=wide;
        by fusion_id;
        var qsim_tester;
        id line;
        run;

/* For each fusion calculate the percent of lines that show bias */
    data flag;
        set CEGS.ase_qsim_tester;
        if qsim_tester ne 0.5 then flag = 1; *when is qsim not equal to 0.5;
        else flag = 0;
        run;

    proc sort data=flag;
        by fusion_id;
        run;

    proc freq data=flag noprint;
        by fusion_id;
        table flag / out=freq;
        run;

    data bias;
        set freq;
        where flag eq 1;
        label percent = ' ';
        rename percent = percent_bias;
        drop count flag;
        run;

/* merge to wide dataset */
    proc sort data=wide;
        by fusion_id;
        run;

    proc sort data=bias;
        by fusion_id;
        run;

    data merged oops;
        merge wide (in=in1) bias (in=in2);
        by fusion_id;
        if in1 and not in2 then percent_bias = 0;
        if in1 then output merged;
        else output oops;
        run;

/* Export for plotting in pythyon */
    proc export data=merged outfile='!MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/qsim_bias_wide.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Summarize to Fusion Level */
proc sort data=flag;
    by fusion_id;
    run;

proc means data=flag;
    by fusion_id;
    output out=sums sum(flag)=num_line_bias;
    run;
    
ods rtf file='!MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/num_lines_bias_freqs.rtf';
proc freq data=sums;
    table num_line_bias;
    run;
ods rtf close;

/* Clean up */
    proc datasets ;
        delete bias;
        delete flag;
        delete freq;
        delete merged;
        delete oops;
        delete wide;
        delete sums;
        run; quit;
