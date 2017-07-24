ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For all events I am going to redo detection with a more-stringent cut-off
   this i will use for later differential detection/gene AS comparisons, etc. */


* Junctions;
data flag_splicing;
  set event.mm10_refseq_splicing_counts;
  where sample_id in ('NSC1','NSC2');
  if apn ge 5 then flag_event_apn_ge5=1;
  else flag_event_apn_ge5=0;

  if apn ge 10 then flag_event_apn_ge10=1;
  else flag_event_apn_ge10=0;

  keep sample_id event_id apn flag_event_apn_ge5 flag_event_apn_ge10;
run;


proc sort data=flag_splicing;
  by event_id;
proc means data=flag_splicing noprint;
  by event_id;
  var flag_event_apn_ge5 flag_event_apn_ge10;
  output out=mean_on mean(flag_event_apn_ge5)=mean_ge5 mean(flag_event_apn_ge10)=mean_ge10;
run;

data flag_event_on;
  set mean_on;
  if mean_ge5 = 1 then do;
     flag_feature_on_ge5=1; flag_feature_ambig_dtct_ge5=0; end;
  else if mean_ge5 = 0 then do;
     flag_feature_on_ge5=0; flag_feature_ambig_dtct_ge5=0; end;
  else do;
     flag_feature_on_ge5=.; flag_feature_ambig_dtct_ge5=1; end;

  if mean_ge10 = 1 then do;
     flag_feature_on_ge10=1; flag_feature_ambig_dtct_ge10=0; end;
  else if mean_ge10 = 0 then do;
     flag_feature_on_ge10=0; flag_feature_ambig_dtct_ge10=0; end;
  else do;
     flag_feature_on_ge10=.; flag_feature_ambig_dtct_ge10=1; end;


  keep event_id flag_feature_on_ge5 flag_feature_ambig_dtct_ge5
  flag_feature_on_ge10 flag_feature_ambig_dtct_ge10;
  rename event_id=feature_id;
run;

* Fragments;
data flag_fragment;
  set event.mm10_refseq_fragment_counts;
  if apn ge 5 then flag_fragment_apn_ge5=1;
  else flag_fragment_apn_ge5=0;

  if apn ge 10 then flag_fragment_apn_ge10=1;
  else flag_fragment_apn_ge10=0;

  keep sample_id fragment_id apn flag_fragment_apn_ge5 flag_fragment_apn_ge10;
run;

proc sort data=flag_fragment;
  by fragment_id;
proc means data=flag_fragment noprint;
  by fragment_id;
  var flag_fragment_apn_ge5 flag_fragment_apn_ge10;
  output out=mean_on mean(flag_fragment_apn_ge5)=mean_ge5 mean(flag_fragment_apn_ge10)=mean_ge10;
run;

data flag_fragment_on;
  set mean_on;
  if mean_ge5 = 1 then do;
     flag_feature_on_ge5=1; flag_feature_ambig_dtct_ge5=0; end;
  else if mean_ge5 = 0 then do;
     flag_feature_on_ge5=0; flag_feature_ambig_dtct_ge5=0; end;
  else do;
     flag_feature_on_ge5=.; flag_feature_ambig_dtct_ge5=1; end;

  if mean_ge10 = 1 then do;
     flag_feature_on_ge10=1; flag_feature_ambig_dtct_ge10=0; end;
  else if mean_ge10 = 0 then do;
     flag_feature_on_ge10=0; flag_feature_ambig_dtct_ge10=0; end;
  else do;
     flag_feature_on_ge10=.; flag_feature_ambig_dtct_ge10=1; end;


  keep fragment_id flag_feature_on_ge5 flag_feature_ambig_dtct_ge5
  flag_feature_on_ge10 flag_feature_ambig_dtct_ge10;
  rename fragment_id=feature_id;
run;

data event.flag_feature_on_all_reps_ge5;
   set flag_event_on flag_fragment_on;
run;


proc freq data=event.flag_feature_on_all_reps_ge5;
  tables flag_feature_on_ge5 flag_feature_on_ge10;
run;

/*

     flag_feature_                             Cumulative    Cumulative
            on_ge5    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0     3183261       96.40       3183261        96.40
                 1      118943        3.60       3302204       100.00

                        Frequency Missing = 56486


     flag_feature_                             Cumulative    Cumulative
           on_ge10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0     3225677       97.36       3225677        97.36
                 1       87413        2.64       3313090       100.00

                        Frequency Missing = 45600

*/

