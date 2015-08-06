/*******************************************************************************
* Filename: check_cegsV_expression_patterns.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Need to check is the number of replicates is driving the
* patterns we are seeing
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname norm '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab virgins from CEGS normalized and centered data */
    data virgin;
        set NORM.ccfus_norm_centered;
        where mating_status eq 'V';
        keep fusion_id line rep uq_log_uq_center;
        run;

    proc sort data=virgin;
        by fusion_id line ;
        run;

/* Average across replicates */
    proc means data=virgin noprint;
        by fusion_id line ;
        output out=means mean(uq_log_uq_center)=mean_exp;
        run;

    data reps;
        set means;
        rename _freq_ = num_reps;
        keep line fusion_id _freq_;
        run;

    * Make permanent dataset;
    data rep_by_line;
        set reps;
        keep line num_reps;
        run;

    proc sort data=rep_by_line nodupkey;
        by line;
        run;

    data SEM.cegsV_line_rep_number;
        set rep_by_line;
        run;

/* Merge rep count with CEGS dataset */
    data cegs;
        set SEM.cegs_virgin_norm_cent;
        keep fusion_id line mean_exp symbol_cat;
        run;

    proc sort data=cegs;
        by line fusion_id;
        run;

    proc sort data=reps;
        by line fusion_id;
        run;

    data merged;
        merge cegs (in=in1) reps (in=in2);
        by line fusion_id;
        run;

/* Create box plots to see if expression is driven by replicate number */
    %macro build_box(gene);
        data curr;
            set merged;
            where symbol_cat ? "&gene";
            run;

        proc sort data=curr;
            by fusion_id num_reps;
            run;

        title "&gene";
        proc sgpanel data=curr;
            panelby fusion_id / onepanel;
            vbox mean_exp / category=num_reps;
            run;
    %mend;

    ods listing close;
    ods pdf file='/home/jfear/mclab/cegs_sem_sd_paper/reports/mean_exp_by_rep.pdf';
    %build_box(Sxl);
    %build_box(dsx);
    %build_box(InR);
    ods pdf close;
    ods listing;

/* Clean Up */
    proc datasets ;
        run; quit;
