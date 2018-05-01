ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* I want to know how the distribution of transcripts by whether they have unique event and by proportion
   of events changes when I use a different APN threshold, AND how the distribution of PB transcripts within
   these groups also changes.

   To do this, I need to redo my APN detection flags for splicing events and for fragments, then group
   transcripts by unique/no unique and into 0-25%,25-50%,50-75%,75-100%,100% bins.

   I should not need to change the event-to-transcript assignment, so I will be using the existing
   feature2xs2gene_exp_only_nomulti dataset.

   I am going to do this is a giant macro so that I can automate the entire process to make this quicker

   Will be using APN>0, APN>=5

 */


%macro binXS(apn,test);
/* 1. Re-flag detection of junctions and fragments */


* Junctions;
data flag_splicing;
  set event.mm10_refseq_splicing_counts_sim;
  where sample_id = "&test.";
  %if &apn.=0 %then %do;
  if apn > &apn. then flag_event_apn_gt0=1;
  else flag_event_apn_gt0=0;
  %end;
  %else %do;
  if apn ge &apn. then flag_event_apn_gt0=1;
  else flag_event_apn_gt0=0;
  %end;
  keep sample_id event_id apn flag_event_apn_gt0;
run;


proc sort data=flag_splicing;
  by event_id;
proc means data=flag_splicing noprint;
  by event_id;
  var flag_event_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data flag_event_on;
  set mean_on;
  length feature_type $8.;
  if mean_gt0 ge 0.5 then flag_feature_on=1;
  else flag_feature_on=0;
  feature_type="junction";
  keep event_id flag_feature_on feature_type;
  rename event_id=feature_id;
run;

* Fragments;
data flag_fragment;
  set event.mm10_refseq_fragment_counts_sim;
  where sample_id ? "&test.";
  %if &apn.=0 %then %do;
  if apn > &apn. then flag_fragment_apn_gt0=1;
  else flag_fragment_apn_gt0=0;
  %end;
  %else %do;
  if apn ge &apn. then flag_fragment_apn_gt0=1;
  else flag_fragment_apn_gt0=0;
  %end;
  keep sample_id fragment_id apn flag_fragment_apn_gt0;
run;

proc sort data=flag_fragment;
  by fragment_id;
proc means data=flag_fragment noprint;
  by fragment_id;
  var flag_fragment_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data flag_fragment_on;
  set mean_on;
  length feature_type $8.;
  if mean_gt0 ge 0.5 then flag_feature_on=1;
  else flag_feature_on=0;
  feature_type="fragment";
  keep fragment_id flag_feature_on feature_type;
  rename fragment_id=feature_id;
run;


data flag_feature_on;
   set flag_event_on flag_fragment_on;
run;

/* 2. For the set of expressed genes without multigene components,
      calculate the proportion of detected events per transcript */

* merge detection flags with transcript_id;
data feat2xs;
  set event.feature2xs2gene_exp_nomult_&test.;
  keep feature_id transcript_id;
run;

data feat_to_keep;
   set event.flagged_feature_short;
   where flag_feature_short=0;
   keep feature_id;
run;

proc sort data=feat2xs nodup;
   by feature_id transcript_id;
proc sort data=feat_to_keep nodup;
   by feature_id;
proc sort data=flag_feature_on nodup;
   by feature_id;
run;

data feat2xs_flag_on;
  merge feat2xs (in=in1) flag_feature_on (in=in2) feat_to_keep (in=in3);
  by feature_id;
  if in1 and in2 and in3;
  flag_feature=1;
run;

* for each transcript, calculate the proportion of events detected ;

proc sort data=feat2xs_flag_on;
   by transcript_id feature_type;
run;

proc means data=feat2xs_flag_on noprint;
  by transcript_id feature_type;
  var flag_feature flag_feature_on;
  output out=num_feat_dtct sum(flag_feature)=num_total_features
                           sum(flag_feature_on)=num_features_dtct;
run;

data junc frag;
   set num_feat_dtct;
   if feature_type="junction" then output junc;
   if feature_type="fragment" then output frag;
   drop _TYPE_ _FREQ_ feature_type;
run;

data junc2;
  set junc;
  rename num_total_features=num_total_junctions num_features_dtct=num_junctions_dtct;
run;

data frag2;
  set frag;
  rename num_total_features=num_total_fragments num_features_dtct=num_fragments_dtct;
run;

proc sort data=junc2;
   by transcript_id;
proc sort data=frag2;
   by transcript_id;
run;

data calc_prop_dtct;
  merge frag2 (in=in1) junc2 (in=in2);
  by transcript_id;
  if not in1 then do;
         num_total_fragments=0;
         num_fragments_dtct=0;
         end;
  if not in2 then do;
         num_total_junctions=0;
         num_junctions_dtct=0;
         end;
  num_total_features=num_total_junctions+num_total_fragments;
  num_features_dtct=num_junctions_dtct+num_fragments_dtct;
  perc_features_dtct=num_features_dtct/num_total_features;
  if num_total_junctions > 0 then do;
  perc_junctions_dtct=num_junctions_dtct/num_total_junctions; end;
  if num_total_fragments > 0 then do;
  perc_fragments_dtct=num_fragments_dtct/num_total_fragments; end;
run;

* flag if transcript has unique event;

data xs_w_uniq;
  set event.flag_xscripts_w_unique;
  keep transcript_id flag_xscript_has_unique;
run;

proc sort data=calc_prop_dtct;
  by transcript_id;
proc sort data=xs_w_uniq;
  by transcript_id;
run;

data xs_prop_dtct_flag_uniq;
  merge calc_prop_dtct (in=in1) xs_w_uniq (in=in2);
  by transcript_id;
  if in1 and in2 then output;
  else if in1 then do;
    flag_xscript_has_unique=0;
    output; end;
run;

/* 3. Bin transcripts by unique/no unique and by proportion events detected */

data bin_xs;
   set xs_prop_dtct_flag_uniq;
   length bin_xscript_perc_dtct $10.;
   if perc_features_dtct = 1 then bin_xscript_perc_dtct = "100%";
   else if perc_features_dtct >= 0.75 then bin_xscript_perc_dtct = "75-99%";
   else if perc_features_dtct >= 0.5 then bin_xscript_perc_dtct = "50-74%";
   else if perc_features_dtct >= 0.25 then bin_xscript_perc_dtct = "25-50%";
   else if perc_features_dtct > 0 then bin_xscript_perc_dtct = "1-24%";
   else bin_xscript_perc_dtct = "0%";
run;

/* 4. Count transcripts per bin */

proc freq data=bin_xs ;
   tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;

/* Make permenant */
data event.bin_xs_by_dtct_apn&apn._&test.;
   set bin_xs;
run;

%mend;


%binXS(0,test1);%binXS(0,test2);
%binXS(1,test1);%binXS(1,test2);
%binXS(2,test1);%binXS(2,test2);
%binXS(5,test1);%binXS(5,test2);
%binXS(10,test1);%binXS(10,test2);
%binXS(25,test1);%binXS(25,test2);


%binXS(0,sim1_test1);%binXS(0,sim1_test2);
%binXS(0,sim2_test1);%binXS(0,sim2_test2);
%binXS(0,sim3_test1);%binXS(0,sim3_test2);

