/* Import unique junction counts, map to distinct junction coordinates, flag on/off */

libname eventloc '/mnt/store/event_sandbox/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import junc counts */

%macro importJunc(sample);

    data WORK.&SAMPLE._JUNC    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/mnt/store/event_sandbox/unique_junction_counts/cvrg_cnts_&SAMPLE..csv" delimiter =
',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $4. ;
       informat seq_name $16. ;
       informat mapped_reads best32. ;
       informat read_length best32. ;
       informat region_length best32. ;
       informat region_depth best32. ;
       informat reads_in_region best32. ;
       informat apn best32. ;
       informat rpkm best32. ;
       informat mean best32. ;
       informat std best32. ;
       informat cv best32. ;
       format sample_id $4. ;
       format seq_name $16. ;
       format mapped_reads best12. ;
       format read_length best12. ;
       format region_length best12. ;
       format region_depth best12. ;
       format reads_in_region best12. ;
       format apn best12. ;
       format rpkm best12. ;
       format mean best12. ;
       format std best12. ;
       format cv best12. ;
    input
              sample_id $
              seq_name  $
              mapped_reads
              read_length
              region_length
              region_depth
              reads_in_region
              apn
              rpkm
              mean
              std
              cv
  ;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  run;

%mend;
%importJunc(NSC1);
%importJunc(NSC2);

/* Merge */

data NSC1_junc2;
   set NSC1_junc;
   keep seq_name apn;
   rename apn=NSC1_apn;
run;


data NSC2_junc2;
   set NSC2_junc;
   keep seq_name apn;
   rename apn=NSC2_apn;
run;

proc sort data=NSC1_junc2;
   by seq_name;
proc sort data=NSC2_junc2;
   by seq_name;
run;

data NSC_junc;
  merge NSC1_junc2 (in=in1) NSC2_junc2 (in=in2);
  by seq_name;
  if in1 and in2;
run;

/* Flag on/off */

data flag_junc_on;
  set NSC_junc;
  if NSC1_apn > 0 or NSC2_apn > 0 then flag_junction_nsc_on=1;
  else flag_junction_nsc_on=0;
run;

proc freq data=flag_junc_on;
  tables flag_junction_nsc_on;
run;

/*
flag_junction_                             Cumulative    Cumulative
        nsc_on    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0     2744894       93.77       2744894        93.77
             1      182421        6.23       2927315       100.00

*/

/* Match to distinct junc coord */

data uniq2coord;
  set evspl.mm10_refseq_junc2uniq_seq;
  drop seq2;
  rename unique_junc_id=seq_name;
run;

proc sort data=uniq2coord;
   by seq_name;
proc sort data=flag_junc_on;
   by seq_name;
run;

data flag_junc_on_coord;
  merge uniq2coord (in=in1) flag_junc_on (in=in2);
  by seq_name;
  if in1 and in2;
run;

data event2coord;
   set evspl.splicing_events_annot_refseq;
   length junction_id $27.;
   junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
   keep junction_id event_id;
run;

proc sort data=event2coord;
   by junction_id;
proc sort data=flag_junc_on_coord;
   by junction_id;
run;

data event2uniq;
  merge event2coord (in=in1) flag_junc_on_coord (in=in2);
  by junction_id;
  if in1 and in2;
run;

/* Get "on" flags for NSCs */

data event_on;
  set event.flag_splicing_on;
  keep flag_event_nsc_on event_id;
run;

proc sort data=event2uniq;
  by event_id;
proc sort data=event_on;
  by event_id;
run;

data event2uniq_on;
  merge event2uniq (in=in1) event_on (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=event2uniq_on;
  tables flag_junction_nsc_on*flag_event_nsc_on;
run;

/*
  flag_junction_nsc_on
            flag_event_nsc_on

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |2811478 |      1 |2811479
           |  93.76 |   0.00 |  93.76
           | 100.00 |   0.00 |
           |  99.74 |   0.00 |
  ---------+--------+--------+
         1 |   7323 | 179652 | 186975
           |   0.24 |   5.99 |   6.24
           |   3.92 |  96.08 |
           |   0.26 | 100.00 |
  ---------+--------+--------+
  Total     2818801   179653  2998454
              94.01     5.99   100.00

Lose one event, but it is an IR event so this shouldn't affect the junction tables
*/

/* Make permenant */

data eventloc.unique_junction2event_mm10;
   set event2uniq_on;
run;


