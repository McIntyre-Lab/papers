/*******************************************************************************
* Filename: qsim_bayesian_prep_data_for_machine.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries

*/

/* Only pull in line*fusions that can be analyzed */
data filter;
    set CEGS.emp_bayesian_input;
    where flag_analyze=1;
    run;

/* Merge with qSIM */
proc sort data=filter;
    by line fusion_id;
    run;

proc sort data=CEGS.ase_qsim_tester;
    by line fusion_id;
    run;

data merged;
    merge filter (in=in1) CEGS.ase_qsim_tester (in=in2);
    by line fusion_id;
    if in1;
    drop apn_: mean_apn sum_tester sum_line sum_both both_:;
    run;

/* Export dataset */
    proc export data=merged outfile="!MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/input/ase_dataset_for_bayesian_w_qsim.csv" dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
    proc datasets nolist;
        delete filter;
        delete merged;
        run; quit;

