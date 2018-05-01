ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Export MISO summary data so I can plot the distribution of PSIs, delta PSIs
   Bayes factors, etc.

   I need the following:
   event_name
   diffs, bayes_factors
   posterior means
   diff flags and bins
*/

data miso_data_for_plots;
   set event.miso_bin_diffs_and_bayes;
   /* I want to calc the CI ranges here too. Do for only NSC and OLD (merged reps) */
   NSC_ci_range=abs(NSC_ci_high-NSC_ci_low);
   OLD_ci_range=abs(OLD_ci_high-OLD_ci_low);
   /* Flag if range is >0.2 */
   if NSC_ci_range > 0.2 then flag_NSC_ci_range_gt02=1;
   else flag_NSC_ci_range_gt02=0;
   if OLD_ci_range > 0.2 then flag_OLD_ci_range_gt02=1;
   else flag_OLD_ci_range_gt02=0;
   drop NSC_ci_low NSC_ci_high OLD_ci_low OLD_ci_high
    NSC1_ci_low NSC1_ci_high NSC2_ci_low NSC2_ci_high 
    OLD1_ci_low OLD1_ci_high OLD2_ci_low OLD2_ci_high 
    NSC_counts NSC_assigned_counts OLD_counts OLD_assigned_counts
    NSC1_counts NSC1_assigned_counts OLD1_counts OLD1_assigned_counts
    NSC2_counts NSC2_assigned_counts OLD2_counts OLD2_assigned_counts;
run;


proc export data=miso_data_for_plots
     outfile="!MCLAB/event_analysis/analysis_output/miso_data_for_plots.csv"
     dbms=csv replace;
run;

