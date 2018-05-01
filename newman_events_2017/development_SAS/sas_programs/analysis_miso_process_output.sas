ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Check: I want to make sure that the posterior means and CIs are the same between comparisons for the same sample
   e.g. NSC1 posterior should be the same regardless of comparison */

data check_miso;
  set event.miso_all_results_nsc_v_old;
    if NSC1_posterior_mean_1 ne . and NSC1_posterior_mean_2 ne . and NSC1_posterior_mean_3 ne . then do;
    if NSC1_posterior_mean_1=NSC1_posterior_mean_2 then flag_NSC1_1v2_okay=1; else flag_NSC1_1v2_okay=0;
    if NSC1_posterior_mean_1=NSC1_posterior_mean_3 then flag_NSC1_1v3_okay=1; else flag_NSC1_1v3_okay=0;
    if NSC1_posterior_mean_2=NSC1_posterior_mean_3 then flag_NSC1_2v3_okay=1; else flag_NSC1_2v3_okay=0;
    end;

    if NSC2_posterior_mean_1 ne . and NSC2_posterior_mean_4 ne . and NSC2_posterior_mean_5 ne . then do;
    if NSC2_posterior_mean_1=NSC2_posterior_mean_4 then flag_NSC2_1v4_okay=1; else flag_NSC2_1v4_okay=0;
    if NSC2_posterior_mean_1=NSC2_posterior_mean_5 then flag_NSC2_1v5_okay=1; else flag_NSC2_1v5_okay=0;
    if NSC2_posterior_mean_4=NSC2_posterior_mean_5 then flag_NSC2_4v5_okay=1; else flag_NSC2_4v5_okay=0;
    end;

    if OLD1_posterior_mean_2 ne . and OLD1_posterior_mean_4 ne . and OLD1_posterior_mean_6 ne . then do;
    if OLD1_posterior_mean_2=OLD1_posterior_mean_4 then flag_OLD1_2v4_okay=1; else flag_OLD1_2v4_okay=0;
    if OLD1_posterior_mean_2=OLD1_posterior_mean_6 then flag_OLD1_2v6_okay=1; else flag_OLD1_2v6_okay=0;
    if OLD1_posterior_mean_4=OLD1_posterior_mean_6 then flag_OLD1_4v6_okay=1; else flag_OLD1_4v6_okay=0;
    end;

    if OLD2_posterior_mean_3 ne . and OLD2_posterior_mean_5 ne . and OLD2_posterior_mean_6 ne . then do;
    if OLD2_posterior_mean_3=OLD2_posterior_mean_5 then flag_OLD2_3v5_okay=1; else flag_OLD2_3v5_okay=0;
    if OLD2_posterior_mean_3=OLD2_posterior_mean_6 then flag_OLD2_3v6_okay=1; else flag_OLD2_3v6_okay=0;
    if OLD2_posterior_mean_5=OLD2_posterior_mean_6 then flag_OLD2_5v6_okay=1; else flag_OLD2_5v6_okay=0;
    end;
run;

proc freq data=check_miso;
   tables flag_NSC1_1v2_okay
          flag_NSC1_1v3_okay
          flag_NSC1_2v3_okay

          flag_NSC2_1v4_okay
          flag_NSC2_1v5_okay
          flag_NSC2_4v5_okay

          flag_OLD1_2v4_okay
          flag_OLD1_2v6_okay
          flag_OLD1_4v6_okay

          flag_OLD2_3v5_okay
          flag_OLD2_3v6_okay
          flag_OLD2_5v6_okay;
run;

/* all okay! */

/* I want to have one set of each posterior for each sample, plus each diff */

data miso_compare_summary;
  set event.miso_all_results_nsc_v_old;
   /* NSC1 */
   if NSC1_posterior_mean_1 ne . then do;
        NSC1_posterior_mean=NSC1_posterior_mean_1;
	NSC1_ci_low=NSC1_ci_low_1;
	NSC1_ci_high=NSC1_ci_high_1;
        NSC1_counts=NSC1_counts_1;
	NSC1_assigned_counts=NSC1_assigned_counts_1;
        end;

   else if NSC1_posterior_mean_2 ne . then do;
        NSC1_posterior_mean=NSC1_posterior_mean_2;
	NSC1_ci_low=NSC1_ci_low_2;
	NSC1_ci_high=NSC1_ci_high_2;
        NSC1_counts=NSC1_counts_2;
	NSC1_assigned_counts=NSC1_assigned_counts_2;
        end;
   else do;
        NSC1_posterior_mean=NSC1_posterior_mean_3;
	NSC1_ci_low=NSC1_ci_low_3;
	NSC1_ci_high=NSC1_ci_high_3;
        NSC1_counts=NSC1_counts_3;
	NSC1_assigned_counts=NSC1_assigned_counts_3;
        end;

   /* NSC2 */
   if NSC2_posterior_mean_1 ne . then do;
        NSC2_posterior_mean=NSC2_posterior_mean_1;
	NSC2_ci_low=NSC2_ci_low_1;
	NSC2_ci_high=NSC2_ci_high_1;
        NSC2_counts=NSC2_counts_1;
	NSC2_assigned_counts=NSC2_assigned_counts_1;
        end;
   else if NSC2_posterior_mean_4 ne . then do;
        NSC2_posterior_mean=NSC2_posterior_mean_4;
	NSC2_ci_low=NSC2_ci_low_4;
	NSC2_ci_high=NSC2_ci_high_4;
        NSC2_counts=NSC2_counts_4;
	NSC2_assigned_counts=NSC2_assigned_counts_4;
        end;
   else do;
        NSC2_posterior_mean=NSC2_posterior_mean_5;
	NSC2_ci_low=NSC2_ci_low_5;
	NSC2_ci_high=NSC2_ci_high_5;
        NSC2_counts=NSC2_counts_5;
	NSC2_assigned_counts=NSC2_assigned_counts_5;
        end;

   /* OLD1 */
   if OLD1_posterior_mean_2 ne . then do;
        OLD1_posterior_mean=OLD1_posterior_mean_2;
	OLD1_ci_low=OLD1_ci_low_2;
	OLD1_ci_high=OLD1_ci_high_2;
        OLD1_counts=OLD1_counts_2;
	OLD1_assigned_counts=OLD1_assigned_counts_2;
        end;
   else if OLD1_posterior_mean_4 ne . then do;
        OLD1_posterior_mean=OLD1_posterior_mean_4;
	OLD1_ci_low=OLD1_ci_low_4;
	OLD1_ci_high=OLD1_ci_high_4;
        OLD1_counts=OLD1_counts_4;
	OLD1_assigned_counts=OLD1_assigned_counts_4;
        end;
   else do;
        OLD1_posterior_mean=OLD1_posterior_mean_6;
	OLD1_ci_low=OLD1_ci_low_6;
	OLD1_ci_high=OLD1_ci_high_6;
        OLD1_counts=OLD1_counts_6;
	OLD1_assigned_counts=OLD1_assigned_counts_6;
        end;

   /* OLD2 */
   if OLD2_posterior_mean_3 ne . then do;
        OLD2_posterior_mean=OLD2_posterior_mean_3;
	OLD2_ci_low=OLD2_ci_low_3;
	OLD2_ci_high=OLD2_ci_high_3;
        OLD2_counts=OLD2_counts_3;
	OLD2_assigned_counts=OLD2_assigned_counts_3;
        end;
   else if OLD2_posterior_mean_5 ne . then do;
        OLD2_posterior_mean=OLD2_posterior_mean_5;
	OLD2_ci_low=OLD2_ci_low_5;
	OLD2_ci_high=OLD2_ci_high_5;
        OLD2_counts=OLD2_counts_5;
	OLD2_assigned_counts=OLD2_assigned_counts_5;
        end;
   else do;
        OLD2_posterior_mean=OLD2_posterior_mean_6;
	OLD2_ci_low=OLD2_ci_low_6;
	OLD2_ci_high=OLD2_ci_high_6;
        OLD2_counts=OLD2_counts_6;
	OLD2_assigned_counts=OLD2_assigned_counts_6;
        end;

   keep event_name NSC1_posterior_mean NSC1_ci_low NSC1_ci_high NSC1_counts NSC1_assigned_counts
        NSC2_posterior_mean NSC2_ci_low	NSC2_ci_high NSC2_counts NSC2_assigned_counts
        OLD1_posterior_mean OLD1_ci_low	OLD1_ci_high OLD1_counts OLD1_assigned_counts
        OLD2_posterior_mean OLD2_ci_low	OLD2_ci_high OLD2_counts
        NSC1_NSC2_diff_1 NSC1_OLD1_diff_2 NSC1_OLD2_diff_3
          NSC2_OLD1_diff_4 NSC2_OLD2_diff_5 OLD1_OLD2_diff_6
          NSC1_NSC2_bayes_factor_1 NSC1_OLD1_bayes_factor_2 NSC1_OLD2_bayes_factor_3
          NSC2_OLD1_bayes_factor_4 NSC2_OLD2_bayes_factor_5 OLD1_OLD2_bayes_factor_6
          NSC_posterior_mean_7 NSC_ci_low_7 NSC_ci_high_7
          OLD_posterior_mean_7 OLD_ci_low_7 OLD_ci_high_7
          NSC_OLD_diff_7 NSC_counts_7 NSC_assigned_counts_7
           OLD_counts_7 OLD_assigned_counts_7 NSC_OLD_bayes_factor_7;
    
      rename NSC1_NSC2_diff_1=NSC1_NSC2_diff
          NSC1_OLD1_diff_2=NSC1_OLD1_diff
          NSC1_OLD2_diff_3=NSC1_OLD2_diff
          NSC2_OLD1_diff_4=NSC2_OLD1_diff
          NSC2_OLD2_diff_5=NSC2_OLD2_diff
          OLD1_OLD2_diff_6=OLD1_OLD2_diff

          NSC1_NSC2_bayes_factor_1=NSC1_NSC2_bayes_factor
          NSC1_OLD1_bayes_factor_2=NSC1_OLD1_bayes_factor
          NSC1_OLD2_bayes_factor_3=NSC1_OLD2_bayes_factor
          NSC2_OLD1_bayes_factor_4=NSC2_OLD1_bayes_factor
          NSC2_OLD2_bayes_factor_5=NSC2_OLD2_bayes_factor
          OLD1_OLD2_bayes_factor_6=OLD1_OLD2_bayes_factor

	  NSC_posterior_mean_7=NSC_posterior_mean
          NSC_ci_low_7=NSC_ci_low
          NSC_ci_high_7=NSC_ci_high
          OLD_posterior_mean_7=OLD_posterior_mean
          OLD_ci_low_7=OLD_ci_low
          OLD_ci_high_7=OLD_ci_high
          NSC_OLD_diff_7=NSC_OLD_diff
          NSC_counts_7=NSC_counts
          NSC_assigned_counts_7=NSC_assigned_counts
          OLD_counts_7=OLD_counts
          OLD_assigned_counts_7=OLD_assigned_counts
          NSC_OLD_bayes_factor_7=NSC_OLD_bayes_factor;
run;

/* Make permenant */

data event.miso_all_results_nsc_v_old_trim;
   set miso_compare_summary;
run;

