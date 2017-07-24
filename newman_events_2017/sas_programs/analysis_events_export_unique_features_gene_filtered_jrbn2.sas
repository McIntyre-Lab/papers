
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Export counts for making plots to visualize the distribution of data post gene-level filtering
   I want to explore the following:
   Specifity*Feature type*Detection (post filtering!) -- thinking about a stacked barplot here!
   Features binned by #num xscripts
   e.g. 1 xs, 2 xs, 3 xs, 4 xs, 5xs, 6-10xs, 11-20xs, >20xs
*/

/* for each feature, get the type, detection (post-filtering), specificity and number of transcripts */

data feat_specificity;
   set event.flag_features_by_specificity_exp;
   keep feature_id xscripts_per_feature flag_feature_unique flag_feature_common flag_feature_constitutive flag_multigene;
run;

data feat_detection;
  set event.features_w_annotations;
  where num_transcripts > 0;
  keep feature_id flag_feature_on feature_type flag_singleton;
run;

data feat_to_keep;
   set event.flagged_feature_short;
   where flag_feature_short=0;
   keep feature_id;
run;

proc sort data=feat_specificity;
   by feature_id;
proc sort data=feat_detection;
   by feature_id;
proc sort data=feat_to_keep;
  by feature_id;
run;

data feat_data_for_plots short_feat;
  merge feat_to_keep (in=in1) feat_detection (in=in2) feat_specificity (in=in3);
  by feature_id;
  if in1 and in2 and in3 then output feat_data_for_plots;
  else if in2 and in3 then output short_feat;
run;

* Export;
proc export data=feat_data_for_plots outfile="!MCLAB/event_analysis/analysis_output/event_analysis_features_retained_with_info_gene_filtered.csv" dbms=csv replace;
run;

data event.feature_data_for_plots_gene_exp;
   set feat_data_for_plots;
run;

