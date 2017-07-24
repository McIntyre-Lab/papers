ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Flag hits with gaps or mismatch */

data flag_partial;
  set event.blast_junc_pb_hits_with_annot;
  if perc_identity >= 95 then flag_perc_ident_ge95=1; else flag_perc_ident_ge95=0;
  if mismatch = 0 then flag_no_mismatch=1; else flag_no_mismatch=0;
  if gapopen = 0 then flag_no_gap=1; else flag_no_gap=0;
run;

proc freq data=flag_partial noprint;
  tables flaG_perc_ident_ge95*flag_no_mismatch*flag_no_gap / out=check;
proc print data=check;
run;

/*
 flag_perc_    flag_no_    flag_no_
 ident_ge95    mismatch       gap       COUNT    PERCENT

      0            0           0          355     0.2278
      0            0           1          426     0.2734
      0            1           0           18     0.0116
      1            0           0          688     0.4415
      1            0           1         4728     3.0340
      1            1           0          777     0.4986
      1            1           1       148842    95.5132

148842 hits to keep
*/


data event_means;
   set event.mean_apn_events_npc;
   if mean_apn_npc > 0 then flag_mean_apn_gt0=1; else flag_mean_apn_gt0=0;
   if mean_apn_npc >= 5 then flag_mean_apn_ge5=1; else flag_mean_apn_ge5=0;
   if mean_apn_npc >= 10 then flag_mean_apn_ge10=1; else flag_mean_apn_ge10=0;
   rename feature_id=event_id;
run;

data event_on;
   set event.features_w_annotations_nomulti;
   where feature_type="splicing" ;
   keep feature_id flag_feature_on flag_junction_annotated flag_intron_retention;
   rename feature_id=event_id;
run;

proc sort data=flag_partial;
   by event_id;
proc sort data=event_means;
   by event_id;
proc sort data=event_on;
   by event_id;
run;

data blast_hits_check_on;
  merge flag_partial (in=in1) event_on (in=in2) event_means;
  by event_id;
  if in1 then flag_event_has_blast=1; else flag_event_has_blast=0;
  if in2 ;
run;


/* add fragment length */

data junc_length;
   set evspl.splicing_events_annot_refseq;
   feature_length=event_size;
   keep event_id feature_length;
run;

proc sort data=blast_hits_check_on;
  by event_id;
proc sort data=junc_length;
  by event_id;
run;

data blast_hits_w_len;
  merge blast_hits_check_on (in=in1) junc_length (in=in2);
  by event_id;
  if in1 and in2 ;
run;


data flag_length;
   set blast_hits_W_len;
   if flag_event_has_blast=1 then do;
      if length >= 0.9 * feature_length then flag_hit_90perc_len=1; else flag_hit_90perc_len=0;
      if length >= 0.95 * feature_length then flag_hit_95perc_len=1; else flag_hit_95perc_len=0;
      if length >= 0.99 * feature_length then flag_hit_99perc_len=1; else flag_hit_99perc_len=0;
      if length >= 1 * feature_length then flag_hit_100perc_len=1; else flag_hit_100perc_len=0;
   end;
run;

/* count -- this is per hit */

proc freq data=flag_length noprint;
  tables flag_event_has_blast*flag_pb_gene_has_refseq*flag_junction_annotated*flag_intron_retention*
         flag_perc_ident_ge95*flag_no_mismatch*flag_no_gap*
         flag_feature_on*flag_hit_90perc_len*flag_mean_apn_gt0*
         flag_mean_apn_ge5 /
         out=pb_hit_count;
run;

proc export data=pb_hit_count outfile="!MCLAB/event_analysis/analysis_output/blast_junctions_to_pacbio_summary_of_hits.csv"
   dbms=csv replace;
run;

/* Extract only hits that are "good" */

data good_blast ;
  set flag_length;
  if flag_pb_gene_has_refseq=1 and flag_perc_ident_ge95=1 and flag_no_mismatch=1 and flag_no_gap=1
  and flag_feature_on=1 and flag_hit_90perc_len=1 ;
  keep event_id;
run;

data all_juncs;
   set flag_length;
   if flag_feature_on=1;
   keep event_id flag_mean_apn_gt0 flag_mean_apn_ge5 flag_junction_annotated flag_intron_retention;
run;

proc sort data=good_blast nodup;
   by event_id;
proc sort data=all_juncs nodup;
  by event_id;
run;


data blast_summary;
  merge all_juncs (in=in1) good_blast (in=in2);
  by event_id;
  if in2 then flag_event_has_blast=1;
  else flag_event_has_blast=0;
run;

proc freq data=blast_summary noprint;
  tables  flag_mean_apn_gt0*flag_mean_apn_ge5*flag_junction_annotated*
          flag_intron_retention*flag_event_has_blast / out=summarize_flags;
run;

proc print data=summarize_flags;
run;

/*

  flag_      flag_
  mean_      mean_     flag_junction_    flag_intron_    flag_event_
 apn_gt0    apn_ge5       annotated        retention      has_blast     COUNT    PERCENT

    1          0              0                0              0          8072     6.4208
    1          0              0                0              1           185     0.1472
    1          0              0                1              0         25492    20.2775
    1          0              0                1              1           590     0.4693
    1          0              1                0              0         39231    31.2061
    1          0              1                0              1         12569     9.9979
    1          1              0                0              0           131     0.1042
    1          1              0                0              1            24     0.0191
    1          1              0                1              0          1143     0.9092
    1          1              0                1              1           267     0.2124
    1          1              1                0              0          7074     5.6270
    1          1              1                0              1         30938    24.6094

*/

/* Make permanent */

data event.blast_junc_to_pb_summary;
  set blast_summary;
run;


data event.blast_junc_to_pb_all_hits;
  set flag_length;
run;
