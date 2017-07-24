ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* I want to know how the distribution of transcripts by whether they have unique event and by proportion
   of events changes when I use a different APN threshold, AND how the distribution of PB transcripts within
   these groups also changes.

   To do this, I need to redo my APN detection flags for splicing events and for fragments, then group
   transcripts by unique/no unique and into 0-25%,25-50%,50-75%,75-100%,100% bins, merge in
   transcripts with PB match and then count for each.

   I should not need to change the event-to-transcript assignment, so I will be using the existing
   feature2xs2gene_exp_only_nomulti dataset.

   I am going to do this is a giant macro so that I can automate the entire process to make this quicker

   Will be using APN>0, APN>=1, APN>=2, APN>=5, APN>=10 

   Going to do this on only the NPC samples, since this is what we have been using

 */


%macro binXS(apn);

/* 1. Re-flag detection of junctions and fragments */


* Junctions;
data flag_splicing;
  set event.mm10_refseq_splicing_counts;
  where sample_id in ('NSC1','NSC2');
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
  if mean_gt0 ge 0.5 then flag_feature_on=1;
  else flag_feature_on=0;
  keep event_id flag_feature_on;
  rename event_id=feature_id;
run;

* Fragments;
data flag_fragment;
  set event.mm10_refseq_fragment_counts;
  where sample_id in ('NSC1','NSC2');
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
  if mean_gt0 ge 0.5 then flag_feature_on=1;
  else flag_feature_on=0;
  keep fragment_id flag_feature_on;
  rename fragment_id=feature_id;
run;


data flag_feature_on;
   set flag_event_on flag_fragment_on;
run;

/* 2. For the set of expressed genes without multigene components,
      calculate the proportion of detected events per transcript */

* merge detection flags with transcript_id;
data feat2xs;
  set event.feature2xs2gene_exp_only_nomulti;
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
   by transcript_id;
run;

proc means data=feat2xs_flag_on noprint;
  by transcript_id;
  var flag_feature flag_feature_on;
  output out=num_feat_dtct sum(flag_feature)=num_total_features
                           sum(flag_feature_on)=num_features_dtct;
run;

data calc_prop_dtct;
  set num_feat_dtct;
  perc_features_dtct=num_features_dtct/num_total_features;
  drop _TYPE_ _FREQ_;
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
  if in1 and in2;
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

/* 5. Subset transcripts with PacBio matches and count per bin */

* Try this one first;
data xs_w_pb;
   set event.pacbio2refseq_id_nomulti;
   keep transcript_id;
run;

proc sort data=xs_w_pb nodup;
  by transcript_id;
proc sort data=bin_xs;
  by transcript_id;
run;

data xs_by_dtct_w_pb;
  merge bin_xs (in=in1) xs_w_pb (in=in2);
  by transcript_id;
  if in2 then flag_has_pacbio=1; else flag_has_pacbio=0;
  if in1;
run;

proc freq data=xs_by_dtct_w_pb;
  where flag_has_pacbio=1;
  tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;


* Now import second list and recount;

     data WORK.SPLICE_MATCH    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/references/pacbio_isoforms_list_for_import.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat pacbio_id $10. ; informat source $6.;    informat feature_type $10. ;
        informat start best32. ;  informat end best32.;   informat frame best32. ;
        informat strand $1. ;     informat score best32.; informat transcript_id $18. ;
        informat match_type $23.; informat note $28. ;
        format pacbio_id $10. ;   format source $6. ;     format feature_type $10. ;
        format start best12. ;    format end best12. ;    format frame best12. ;
        format strand $1. ;       format score best12. ;  format transcript_id $18. ;
        format match_type $23. ;  format note $28. ;
        input pacbio_id $ source $ feature_type $ start end
              frame strand $ score transcript_id $ match_type $ note $
              ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data pb2refseq;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep pacbio_id splice_match_id match_type;
   rename splice_match_id=transcript_id;
run; *8722;

data pb2keep;
  set event.pacbio2refseq_id_nomulti;
  keep pacbio_id;
run;

proc sort data=pb2keep nodup;
  by pacbio_id;
proc sort data=pb2refseq nodup;
  by pacbio_id transcript_id;
run;

data pb2refseq2;
  merge pb2refseq (in=in1) pb2keep (in=in2);
  by pacbio_id;
  if in1 and in2;
  drop pacbio_id;
run;

data xs_list_1;
   set pb2refseq;
   keep transcript_id;
run;

data xs_list_2;
   set pb2refseq2;
   keep transcript_id;
run;

proc sort data=xs_list_1 nodup;
  by transcript_id;
proc sort data=xs_list_2 nodup;
  by transcript_id;
proc sort data=bin_xs;
  by transcript_id;
run;

data bin_xs_list_1;
  merge bin_xs (in=in1) xs_list_1 (in=in2);
  by transcript_id;
  if in2 then flag_has_pacbio=1; else flag_has_pacbio=0;
  if in1;
run;

proc freq data=bin_xs_list_1;
  where flag_has_pacbio=1;
  tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;

data bin_xs_list_2;
  merge bin_xs (in=in1) xs_list_2 (in=in2);
  by transcript_id;
  if in2 then flag_has_pacbio=1; else flag_has_pacbio=0;
  if in1;
run;

proc freq data=bin_xs_list_2;
  where flag_has_pacbio=1;
  tables flag_xscript_has_unique*bin_xscript_perc_dtct;
run;

/* Make permenant */

data event.bin_xscripts_by_dtct_apn&apn.;
   set bin_xs_list_1;
run;


%mend;


%binXS(0);
%binXS(1);
%binXS(2);
%binXS(5);
%binXS(10);
%binXS(25);

