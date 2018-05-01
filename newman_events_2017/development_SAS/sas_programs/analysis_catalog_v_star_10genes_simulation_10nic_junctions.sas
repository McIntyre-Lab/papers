/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 59-gene simulation, I want to count the number of times that Event and STAR
   detect the 9 NIC junctions that were simulated

   There were 10 NIC junctions, but I messed up one of the NIC transcripts, so I
   am going to exclude this one*/


data nic_sim;
   format junction_id $50.;
   input junction_id $ ;
   datalines;
   chr6:147719923:147727642:-
   chr6:119921643:119922471:+
   chr5:53466228:53600832:+
   chr16:36875263:36885010:+
   chr3:88736088:88748781:-
   chr5:129109661:129128512:+
   chr7:4631022:4631122:-
   chr7:132837207:132850590:-
   chr13:106823051:106836211:-
   ;
run;

/* How many of these junctions are detected with STAR and via the catalog? */

data cat_junc;
   set event.catalog_junctions_10genes;
   if apn>0 then flag_detected=1; else flag_detected=0;
   keep sample_id event_id flag_detected;
   rename event_id=seq_name;
run;

proc sort data=cat_junc;
   by seq_name sample_id;
proc transpose data=cat_junc out=cat_junc_sbys;
   by seq_name;
   id sample_id;
   var flag_detected;
run;

data cat_seq2id;
  set eventloc.unique_junction2event_mm10;
  keep seq_name junction_id;
run;

proc sort data=cat_junc_sbys;
  by seq_name;
proc sort data=cat_seq2id nodup;
  by seq_name;
run;

data cat_junc_id;
  merge cat_junc_sbys (in=in1) cat_seq2id (in=in2);
  by seq_name;
  if in1 and in2;
  drop seq_name;
run;

proc sort data=cat_junc_id;
  by junction_id;
proc sort data=nic_sim;
  by junction_id;
run;

data cat_junc_id_sim;
  merge nic_sim (in=in1) cat_junc_id (in=in2);
  by junction_id;
  if in1 and in2 then output;
  else if in1 then do;
    sample_01=0; sample_02=0; sample_03=0;
    sample_04=0; sample_05=0; sample_06=0;
    output;
    end;
  drop _NAME_;
  rename sample_01=flag_sample1_events_dtct
         sample_02=flag_sample2_events_dtct
         sample_03=flag_sample3_events_dtct
         sample_04=flag_sample4_events_dtct
         sample_05=flag_sample5_events_dtct
         sample_06=flag_sample6_events_dtct;
run;

/* STAR junctions */


data star_junc;
  set event.star_junctions_10genes;
  length junction_id $50.;
  if max_overhang lt 16 then delete;
  if strand=1 then strand_str="+"; else strand_str="-";
  donor_site=intron_start-1;
  acceptor_site=intron_stop;
  junction_id=catx(":",chr,donor_site,acceptor_site,strand_str);
  if num_unique_mapped_reads > 0 then flag_detected=1;
  else flag_detected=0;
  keep sample_id junction_id flag_detected;
run;

proc sort data=star_junc;
   by junction_id sample_id;
proc transpose data=star_junc out=star_junc_sbys;
   by junction_id;
   id sample_id;
   var flag_detected;
run;

data star_junc_sbys2;
  set star_junc_sbys;
  drop _NAME_;
  rename sample_01=flag_sample1_star_dtct
         sample_02=flag_sample2_star_dtct
         sample_03=flag_sample3_star_dtct
         sample_04=flag_sample4_star_dtct
         sample_05=flag_sample5_star_dtct
         sample_06=flag_sample6_star_dtct;
run;

data star_junc_sbys3;
  set star_junc_sbys2;
  array change _numeric_;
  do over change;
   if change=. then change=0;
  end;
run;



proc sort data=star_junc_sbys3 nodup;
  by junction_id;
proc sort data=cat_junc_id_sim nodup;
  by junction_id;
run;

data cat2star_sim;
  merge cat_junc_id_sim (in=in1) star_junc_sbys3 (in=in2);
  by junction_id;
  if in1 and in2 then output;
  else if in1 then do;
       flag_sample1_star_dtct=0;
       flag_sample2_star_dtct=0;
       flag_sample3_star_dtct=0;
       flag_sample4_star_dtct=0;
       flag_sample5_star_dtct=0;
       flag_sample6_star_dtct=0;
       output; end;
run;


/* Flag if detected in all samples */

data cat2star_sim_flag;
  set cat2star_sim;
  if sum(flag_sample1_events_dtct,flag_sample2_events_dtct,flag_sample3_events_dtct,
         flag_sample4_events_dtct,flag_sample5_events_dtct,flag_sample6_events_dtct) = 6
     then flag_all_events_dtct=1; else flag_all_events_dtct=0;

  if sum(flag_sample1_star_dtct,flag_sample2_star_dtct,flag_sample3_star_dtct,
         flag_sample4_star_dtct,flag_sample5_star_dtct,flag_sample6_star_dtct) = 6
     then flag_all_star_dtct=1; else flag_all_star_dtct=0;

  if sum(flag_sample1_events_dtct,flag_sample2_events_dtct,flag_sample3_events_dtct,
         flag_sample4_events_dtct,flag_sample5_events_dtct,flag_sample6_events_dtct) > 0
     then flag_any_events_dtct=1; else flag_any_events_dtct=0;

  if sum(flag_sample1_star_dtct,flag_sample2_star_dtct,flag_sample3_star_dtct,
         flag_sample4_star_dtct,flag_sample5_star_dtct,flag_sample6_star_dtct) > 0
     then flag_any_star_dtct=1; else flag_any_star_dtct=0;
run;

/* Count */

proc freq data=cat2star_sim_flag noprint;
  tables flag_sample1_events_dtct*flag_sample1_star_dtct / out=sample1_events_v_star;
  tables flag_sample2_events_dtct*flag_sample2_star_dtct / out=sample2_events_v_star;
  tables flag_sample3_events_dtct*flag_sample3_star_dtct / out=sample3_events_v_star;
  tables flag_sample4_events_dtct*flag_sample4_star_dtct / out=sample4_events_v_star;
  tables flag_sample5_events_dtct*flag_sample5_star_dtct / out=sample5_events_v_star;
  tables flag_sample6_events_dtct*flag_sample6_star_dtct / out=sample6_events_v_star;
  tables flag_all_events_dtct*flag_all_star_dtct / out=all_events_v_star;
  tables flag_any_events_dtct*flag_any_star_dtct / out=any_events_v_star;
run;

proc print data=sample1_events_v_star;
run;
proc print data=sample2_events_v_star;
run;
proc print data=sample3_events_v_star;
run;
proc print data=sample4_events_v_star;
run;
proc print data=sample5_events_v_star;
run;
proc print data=sample6_events_v_star;
run;
proc print data=all_events_v_star;
run;
proc print data=any_events_v_star;
run;

/* Tables 

SAMPLE 1:
 flag_sample1_    flag_sample1_
  events_dtct       star_dtct      COUNT

       1                0            3
       1                1            6


SAMPLE 2:

flag_sample2_    flag_sample2_
 events_dtct       star_dtct      COUNT

      1                0            3
      1                1            6


SAMPLE 3:
flag_sample3_    flag_sample3_
 events_dtct       star_dtct      COUNT

      1                0            3
      1                1            6

SAMPLE 4:
 flag_sample4_    flag_sample4_
  events_dtct       star_dtct      COUNT    PERCENT

       1                0            3      33.3333
       1                1            6      66.6667

SAMPLE 5:
flag_sample5_    flag_sample5_
 events_dtct       star_dtct      COUNT

      1                0            3
      1                1            6

SAMPLE 6:
flag_sample6_    flag_sample6_
 events_dtct       star_dtct      COUNT

      1                0            3
      1                1            6

ALL SAMPLES:
flag_all_
 events_     flag_all_
   dtct      star_dtct    COUNT

    1            0          3
    1            1          6

ANY SAMPLE:
 flag_any_
  events_     flag_any_
    dtct      star_dtct    COUNT

     1            0          3
     1            1          6


*/

/* Junctions missed by STAR */

proc print data=cat2star_sim_flag;
  where flag_any_star_dtct=0;
run;

/*
chr5:129109661:129128512:+
chr6:119921643:119922471:+
chr6:147719923:147727642:-
*/


