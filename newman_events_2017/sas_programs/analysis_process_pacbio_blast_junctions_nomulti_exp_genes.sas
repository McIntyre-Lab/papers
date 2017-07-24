ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that junctions are BLASTing to the correct gene. For those that don't flag and put into a separate
   dataset of further examination

1. Check that events are detected
2. Remove hits that are too short, flag if 95+%, 99+% or 100% of original feature length
3. Flag if events are going to the correct gene or not

*/

data blast_hits;
  set event.blast_dtct_junc2pb_nomult;
run;

data feature_on;
   set event.features_w_annotations_nomulti;
   keep feature_id flag_feature_on;
   rename feature_id=event_id;
run;

proc sort data=blast_hits;
   by event_id;
proc sort data=feature_on;
   by event_id;
run;

data blast_hits_check_on not_exp_gene;
  merge blast_hits (in=in1) feature_on (in=in2);
  by event_id;
  if in1 and in2 then output blast_hits_check_on;
  else if in1 then output not_exp_gene;
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
  if in1 and in2;
run;


*drop if hit <90% of length, flag if 95, 99, 100 perc;

data flag_length;
   set blast_hits_W_len;
   if length < 0.9 * feature_length then delete;
   if length >= 0.9 * feature_length then flag_hit_90perc_len=1; else flag_hit_90perc_len=0;
   if length >= 0.95 * feature_length then flag_hit_95perc_len=1; else flag_hit_95perc_len=0;
   if length >= 0.99 * feature_length then flag_hit_99perc_len=1; else flag_hit_99perc_len=0;
   if length >= 1 * feature_length then flag_hit_100perc_len=1; else flag_hit_100perc_len=0;
run;

/*
 flag_hit_90perc_                             Cumulative    Cumulative
              len    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                1       96275      100.00         96275       100.00


 flag_hit_95perc_                             Cumulative    Cumulative
              len    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         207        0.22           207         0.22
                1       96068       99.78         96275       100.00


 flag_hit_99perc_                             Cumulative    Cumulative
              len    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         402        0.42           402         0.42
                1       95873       99.58         96275       100.00


        flag_hit_                             Cumulative    Cumulative
      100perc_len    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         402        0.42           402         0.42
                1       95873       99.58         96275       100.00
*/

/* Flag if event is going to the correct gene */

data flag_gene;
  set flag_length;
  if count(refseq_gene_id,gene_id) > 0 then flag_correct_gene=1;
  else flag_correct_gene=0;
run;

proc freq data=flag_gene;
  tables flag_correct_gene flag_correct_gene*flag_junction_annotated;
run;

/*

  flag_correct_                             Cumulative    Cumulative
           gene    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0        3803        3.95          3803         3.95
              1       92472       96.05         96275       100.00


 flag_correct_gene
           flag_junction_annotated

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    490 |   3313 |   3803
          |   0.51 |   3.44 |   3.95
          |  12.88 |  87.12 |
          |  16.13 |   3.55 |
 ---------+--------+--------+
        1 |   2548 |  89924 |  92472
          |   2.65 |  93.40 |  96.05
          |   2.76 |  97.24 |
          |  83.87 |  96.45 |
 ---------+--------+--------+
 Total        3038    93237    96275
              3.16    96.84   100.00


So far so good. I want to collapse these to individual events
Basically, if an event goes to at least its correct gene, then keep. Else, drop

*/

proc sort data=flag_gene;
  by event_id;
proc means data=flag_gene noprint;
  by event_id;
  var flag_pb_gene_has_refseq flag_junction_annotated flag_intron_retention flag_hit_99perc_len
      flag_hit_95perc_len flag_feature_on flag_correct_gene;
   output out=summarize_blast_hit_flags max=;
run;

proc freq data=summarize_blast_hit_flags noprint;
  tables flag_pb_gene_has_refseq*flag_junction_annotated*flag_intron_retention*flag_hit_99perc_len*
         flag_hit_95perc_len*flag_feature_on*flag_correct_gene / out=summarize_flags;
run;

proc print data=summarize_flags;
run;

/*
flag_pb_                                                     flag_    flag_
ene_has_ flag_junction_ flag_intron_  flag_hit_  flag_hit_ feature_ correct_
 refseq     annotated     retention  99perc_len 95perc_len    on      gene   COUNT

GOOD:
   1            0             0           0          0         1        1       14
   1            0             1           0          1         1        1        3
   1            0             0           0          1         1        1        8
   1            0             0           1          1         1        1     1021
   1            0             1           0          0         1        1        4
   1            0             1           1          1         1        1      797
   1            1             0           0          0         1        1       40
   1            1             0           0          1         1        1       45
   1            1             0           1          1         1        1    42402


BAD:


   0            0             0           0          1         1        0        1
   0            0             0           1          1         1        0       50
   0            0             1           0          0         1        0        4
   0            0             1           0          1         1        0        1
   0            0             1           1          1         1        0      181
   0            1             0           0          0         1        0        5
   0            1             0           0          1         1        0        3
   0            1             0           1          1         1        0     1768
   1            0             0           0          1         1        0        1
   1            0             0           1          1         1        0        6
   1            0             1           0          0         1        0        2
   1            0             1           0          1         1        0        3
   1            0             1           1          1         1        0       48
   1            1             0           0          0         1        0        2
   1            1             0           0          1         1        0        3
   1            1             0           1          1         1        0      174

*/

/* Make permenant */

data event.blast_dtct_jnc2pb_nomult_flags;
   set flag_gene;
run;

data event.blast_dtct_jnc2pb_nomult_summary;
   set summarize_blast_hit_flags;
run;


