/*******************************************************************************
* Filename: emp_bayesian_prep_data_for_machine.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname MYCEGS '!HOME/storage/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Removes line*mv that have less than 3 */
    /* Count How Many reps are present */
        * The bayesian script requires 3 replicates so I need to select these;
        proc sort data=all_ase nodupkey out=uniq;
            by line mating_status rep;
            run; * Uniq is a list of samples and reps;

        proc freq data=uniq noprint;
            table line*mating_status/out=counts;
            run; 
        
    /* Only keep if >= 3 reps */
        data has_3_reps;
            set counts;
            if count ge 3;
            keep line mating_status ;
            run; * 136 lines*mv have 3 or more reps.;

        * Merge on to all_ase and only keeps samples with 3 reps;
        proc sort data=has_3_reps;
            by line mating_status;
            run;
            
        proc sort data=all_ase;
            by line mating_status;
            run;

        data merged;
            merge all_ase (in=in1) has_3_reps (in=in2);
            by line mating_status;
            if in2;
            run;

        * I had some issues where there were 3 reps, but name 4 5 6, so use a count instead of rep;
        proc sort data=merged;
            by line mating_status fusion_id;
            run;

        data merged_cnt;
            set merged;
            count + 1;
            if first.fusion_id then count = 1;
            by line mating_status fusion_id;
            run;

        proc sort data=merged_cnt;
            by DESCENDING count;
            run; * r208 V has 6 reps;

    /* Transpose dataset so reps in columns */
        proc sort data=merged_cnt;
            by line mating_status fusion_id;
            run;

        proc transpose data=merged_cnt out=LINE_flip prefix=LINE_TOTAL_;
            by line mating_status fusion_id;
            var LINE_TOTAL;
            id count;
            run;
            
        proc transpose data=merged_cnt out=TESTER_flip prefix=TESTER_TOTAL_;
            by line mating_status fusion_id;
            var TESTER_TOTAL;
            id count;
            run;

        proc transpose data=merged_cnt out=TESTER_flip prefix=BOTH_TOTAL_;
            by line mating_status fusion_id;
            var BOTH_TOTAL;
            id count;
            run;

        proc transpose data=merged_cnt out=BOTH_flip prefix=BOTH_TOTAL_;
            by line mating_status fusion_id;
            var BOTH_TOTAL;
            id count;
            run;

        data wide_merge;
            retain line mating_status fusion_id 
            LINE_TOTAL_1 LINE_TOTAL_2 LINE_TOTAL_3 LINE_TOTAL_4 LINE_TOTAL_5 LINE_TOTAL_6 
            TESTER_TOTAL_1 TESTER_TOTAL_2 TESTER_TOTAL_3 TESTER_TOTAL_4 TESTER_TOTAL_5 TESTER_TOTAL_6
            BOTH_TOTAL_1 BOTH_TOTAL_2 BOTH_TOTAL_3 BOTH_TOTAL_4 BOTH_TOTAL_5 BOTH_TOTAL_6 
            ;
            format LINE_TOTAL_1 LINE_TOTAL_2 LINE_TOTAL_3 LINE_TOTAL_4 LINE_TOTAL_5 LINE_TOTAL_6 
            TESTER_TOTAL_1 TESTER_TOTAL_2 TESTER_TOTAL_3 TESTER_TOTAL_4 TESTER_TOTAL_5 TESTER_TOTAL_6
            BOTH_TOTAL_1 BOTH_TOTAL_2 BOTH_TOTAL_3 BOTH_TOTAL_4 BOTH_TOTAL_5 BOTH_TOTAL_6
            best12.;
            merge LINE_flip (drop=_name_) TESTER_FLIP(drop=_name_) BOTH_FLIP(drop=_name_);
            by line mating_status fusion_id;
            sum_line = sum(LINE_TOTAL_1, LINE_TOTAL_2, LINE_TOTAL_3, LINE_TOTAL_4, LINE_TOTAL_5, LINE_TOTAL_6);
            sum_tester = sum(TESTER_TOTAL_1, TESTER_TOTAL_2, TESTER_TOTAL_3, TESTER_TOTAL_4, TESTER_TOTAL_5, TESTER_TOTAL_6);
            sum_both = sum(BOTH_TOTAL_1, BOTH_TOTAL_2, BOTH_TOTAL_3, BOTH_TOTAL_4, BOTH_TOTAL_5, BOTH_TOTAL_6);
            sum_total = sum(sum_line, sum_tester, sum_both);
            run;
        
/* Create flag_analyze for each fusion if APN > 0 for at least 1 rep */
    * A fusion is called analyzable if it has an APN > 0 for at least one rep;

    * Merge on fusion lengths so that APN can be estimated from total_reads_counted;
    data fusion_length;
        set DMEL551.fb551_si_fusions;
        len = end - start + 1;
        keep fusion_id len;
        run;

    proc sort data=fusion_length nodups;
        by fusion_id;
        run;

    proc sort data=merged;
        by fusion_id;
        run;

    data merge_len;
        merge merged (in=in1) fusion_length (in=in2);
        by fusion_id;
        if in1;
        run;

    * Estimate APN;
    data apn;
        set merge_len;
        APN = (total_reads_counted * 96) / len; * read length is ~96bp;
        if APN > 0 then flag_apn = 1; else flag_apn =0;
        keep line mating_status rep fusion_id flag_apn apn;
        run;

    proc sort data=apn;
        by line mating_status fusion_id;
        run;
        
    proc transpose data=apn out=apnT prefix=apn_;
        by line mating_status fusion_id;
        var apn;
        id rep;
        run;

    proc means data=apn noprint;
        by line mating_status fusion_id;
        output out=apnm mean(apn)=mean_apn;
        run;
        

    * Sum reps to make flag;
    proc means data=apn noprint;
        by line mating_status fusion_id;
        output out=sums sum(flag_apn)=sums;
        run;

    data flag_analyze;
        set sums;
        if sums > 0 then flag_analyze = 1; else flag_analyze = 0;
        keep line mating_status fusion_id flag_analyze;
        run;

    proc freq data=flag_analyze noprint;
        table line*mating_status*flag_analyze /out=myfreq;
        run;

    proc sort data=myfreq;
        by mating_status;
        run;

    proc means data=myfreq;
        by mating_status;
        var count;
        run;
    
/* Package dataset with flags and counts */
    proc sort data=wide_merge;
        by line mating_status fusion_id;
        run;

    proc sort data=flag_analyze;
        by line mating_status fusion_id;
        run;

    proc sort data=apnt;
        by line mating_status fusion_id;
        run;

    proc sort data=apnm;
        by line mating_status fusion_id;
        run;

    data big_data;
        merge wide_merge (in=in1) flag_analyze (in=in2) apnt (in=in3 drop=_name_) apnm (in=in4 drop=_type_ _freq_);
        by line mating_status fusion_id;
        if in1;
        if sum_line = 0 and sum_tester = 0 then flag_analyze = 0;
        run;
    
    data CEGS.emp_bayesian_input;
        set big_data;
        run;

/* Export full dataset */
    data out;
        set CEGS.emp_bayesian_input;
        drop apn_: mean_apn sum_tester sum_line sum_both BOTH_:;
        run;
        
    proc export data=out outfile="!MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/input/ase_dataset_for_bayesian.csv" dbms=csv replace;
    putnames=yes;
    run;

/* Freqs */
    /* 
    proc sort data=CEGS.emp_bayesian_input;
        by line mating_status;
        run;

    proc freq data=CEGS.emp_bayesian_input;
        by line mating_status;
        tables flag_analyze;
        run;

    */

/* Clean up */
    proc datasets;
        delete all_ase;
        delete apn;
        delete apnt;
        delete apnm;
        delete big_data;
        delete counts;
        delete flag_analyze;
        delete fusions;
        delete fusions2;
        delete fusion_length;
        delete has_3_reps;
        delete lines;
        delete line_flip;
        delete both_flip;
        delete merged;
        delete merged_cnt;
        delete merge_len;
        delete myfreq;
        delete parms;
        delete regstry;
        delete sums;
        delete tester_flip;
        delete uniq;
        delete wide_merge;
        delete out;
        run;quit;
