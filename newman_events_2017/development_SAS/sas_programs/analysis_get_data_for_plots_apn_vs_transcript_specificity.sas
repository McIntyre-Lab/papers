ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* I want to take singletons, fragments and junctions and output their transcript-specificity
   and mean APN to make some plots -- do for detected events from detected genes */

data spec_flags;
   set event.flag_features_by_specificity_exp;
   keep feature_id flag_feature_unique flag_feature_common flag_feature_constitutive;
run;

data feature_types;
   set event.features_w_annotations;
   keep feature_id feature_type flag_singleton;
run;

data fragments_on;
   set event.flag_fragment_on;
   where flag_fragment_nsc_on=1;
   keep fragment_id;
   rename fragment_id=feature_id;
run;

data splicing_on;
   set event.flag_splicing_on;
   where flag_event_nsc_on=1;
   keep event_id;
   rename event_id=feature_id;
run;

data features_on;
   set splicing_on fragments_on;
run;


data splicing_counts;
   set event.mm10_refseq_splicing_counts;
   where sample_id ? "NSC";
   keep sample_id event_id apn;
   rename event_id=feature_id;
run;


data fragment_counts;
   set event.mm10_refseq_fragment_counts;
   where sample_id ? "NSC";
   keep sample_id fragment_id apn;
   rename fragment_id=feature_id;
run;

data feature_counts;
   set splicing_counts fragment_counts;
run;

proc sort data=feature_counts;
   by feature_id;
proc means data=feature_counts noprint;
   by feature_id;
   var apn;
   output out=mean_feature_count mean=;
run;

proc sort data=mean_feature_count;
   by feature_id;
proc sort data=features_on;
   by feature_id;
proc sort data=feature_types;
   by feature_id;
proc sort data=spec_flags;
   by feature_id;
run;

data features_all_data;
   merge feature_types (in=in1) features_on (in=in2) spec_flags (in=in3) mean_feature_count (in=in4);
   by feature_id;
   if in2 and in3;
   drop _TYPE_ _FREQ_;
run;

proc export data=features_all_data outfile="!MCLAB/event_analysis/analysis_output/event_analysis_counts_by_feature_type_specificity.csv"
   dbms=csv replace;
run;


