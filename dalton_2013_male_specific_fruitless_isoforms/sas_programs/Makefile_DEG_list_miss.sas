/********************************************************************************
* Create induced and repressed gene list for the different treatment types
* based on ANOVA results.
*
* Michelle had done this on her side, but I created all of these to test and
* make sure there were no errors. I compared my output list with the ones she
* sent and everyhing matched up.
********************************************************************************/

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Create FruMA Male Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2_miss
    * 
    * DATASET: FRU.fruA_male_fus_ind_rep_miss
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruA_male.sas';

/* Look at differences in fusions */
    data miss;
        set FRU.fruA_male_fusions_ind_rep_miss;
        rename flag_fruA_male_ind = miss_ind;
        rename flag_fruA_male_rep = miss_rep;
        run;

    data truncated;
        set FRU.fruA_male_fusions_ind_rep;
        rename flag_fruA_male_ind = trunc_ind;
        rename flag_fruA_male_rep = trunc_rep;
        run;

    proc sort data=miss;
        by fusion_id;
        run;

    proc sort data=truncated;
        by fusion_id;
        run;

    data merged oops;
        merge truncated (in=in1) miss (in=in2);
        by fusion_id;
        run;

/* Look at differences in FBgns */
    data miss;
        set FRU.fruA_male_FBgn_ind_rep_miss;
        rename flag_ind = miss_ind;
        rename flag_rep = miss_rep;
        run;

    data truncated;
        set FRU.fruA_male_FBgn_ind_rep;
        rename flag_ind = trunc_ind;
        rename flag_rep = trunc_rep;
        run;

    proc sort data=miss;
        by FBgn_cat;
        run;

    proc sort data=truncated;
        by FBgn_cat;
        run;

    data both differ ;
        merge truncated (in=in1) miss (in=in2);
        by FBgn_cat;
        if in1 and in2 then output both;
        else output differ;
        run;

    data differ2;
        set differ;
        rename FBgn_cat = primary_FBgn;
        run;

    proc sort data=differ2;
        by primary_fbgn;
        run;

    proc sort data=DMEL530.symbol2fbgn;
        by primary_fbgn;
        run;

    data differ_w_name;
        merge differ2 (in=in1) DMEL530.symbol2fbgn (in=in2);
        by primary_FBgn;
        if in1 and in2;
        run;

    data diff_list;
        set differ_w_name;
        keep symbol;
        run;

    proc sort data=diff_list;
        by symbol;
        run;

    proc export data=diff_list
            outfile= '/home/jfear/mclab/arbeitman/arbeitman_fru_network/reports_internal/genes_changed_w_missing.txt'
            dbms=CSV replace;
            putnames=no;
            run;
