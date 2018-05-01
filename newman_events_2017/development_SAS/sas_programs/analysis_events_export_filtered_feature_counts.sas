ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count number of remaining exons/fragments and export for making plots. Add in number of transcripts per feature so we can look at distribution.

Output should look the same as the program "analysis_events_export_feature_counts.sas", so I can just
pull the final dataset from there and append filtering flags */

* Fusions and fragments to filter;

data short_frag;
 set event.flagged_fragment_length;
 rename fragment_id=feature_id flag_fragment_lt_min_bp=flag_feature_short;
run;

data short_fus;
  set event.flagged_fusion_length;
 rename fusion_id=feature_id flag_fusion_lt_min_bp=flag_feature_short;
run;

data short_event;
  set event.flagged_event_length;
  drop flag_event_lt_min_bp;
*since we align to events directly, I want to drop events that are shorter than the read length;
 rename event_id=feature_id flag_event_short=flag_feature_short;
run;

* stack;

data feature_flagged_short;
  length feature_id $450.;
  format feature_id $450.;
  set short_event short_frag short_fus;
run;


* Get annotated features;
data features_annot;
   set event.features_w_annotations;
run;

proc sort data=features_annot;
  by feature_id;
proc sort data=feature_flagged_short;
  by feature_id;
run;

data features_annot_flag_short;
  merge features_annot (in=in1) feature_flagged_short (in=in2) ;
  by feature_id;
  if in1 and in2;
run;

/* Make permenant the stacked flags -- I can always remerge with annotations */

data event.flagged_feature_short;
   set feature_flagged_short;
run;

/* Export for plots */

proc export data=features_annot_flag_short 
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_features_w_annotations_flagged_short.csv" dbms=csv replace;
run;




