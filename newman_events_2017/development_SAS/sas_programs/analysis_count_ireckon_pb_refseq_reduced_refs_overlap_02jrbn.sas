/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the mouse NPC data, compare iReckon to:
   (1) PacBio
   (2) Refseq
   (3) Reduced reference: 100% events, APN>0
   (4) Reduced reference: 75% events, APN>5

   I am doing this for each NPC sample separately as the transcript from iReckon don't necessarily match
   between replicates

   For iReckon I will add in the rpkm_bin information as well:
     bin 0: log_rpkm=0
     bin 1: log_rpkm=0 -> 0.1
     bin 2: log_rpkm=0.1 -> 0.5
     bin 3: log_rpkm=0.5 -> 1.5
     bin 4: log_rpkm> 1.5
*/

/* iReckon transcripts */

data ireckon;
   set event.ireckon_results_nsc;
   length ir_transcript_id $100.;
   ir_transcript_id=catx("_",gene_id,transcript_id);
   keep ir_transcript_id sample_id;
run;

/* iReckon to Refseq BLAST best hits */

data ir2refseq;
   set event.ireckon2refseq_blast_best;
   keep query_id sample transcript_id;
   rename query_id=ir_transcript_id sample=sample_id;
run;

/*quick check: is there a unique entry fror each ireckon transcript? */
proc sort data=ir2refseq;
   by sample_id ir_transcript_id;
proc freq data=ir2refseq noprint;
  by sample_id ;
  tables ir_transcript_id / out=rs_per_ir;
run;

data rs_per_ir2;
  set rs_per_ir;
  if count > 1;
run; *good: 1 match per ireckon;

proc sort data=ir2refseq;
   by sample_id ir_transcript_id;
proc sort data=ireckon;
   by sample_id ir_transcript_id;
run;

data ireckon_refseq_best;
  merge ireckon (in=in1) ir2refseq (in=in2);
  by sample_id ir_transcript_id;
  if in2 then flag_has_refseq_hit=1; else flag_has_refseq_hit=0;
run;

/* iReckon to PacBio BLAST best hits */


data ir2pacbio;
   set event.ireckon2pacbio_blast_best;
   keep query_id sample transcript_id;
   rename query_id=ir_transcript_id sample=sample_id transcript_id=pacbio_id;
run;

/*quick check: is there a unique entry fror each ireckon transcript? */
proc sort data=ir2pacbio;
   by sample_id ir_transcript_id;
proc freq data=ir2pacbio noprint;
  by sample_id ;
  tables ir_transcript_id / out=pb_per_ir;
run;

data pb_per_ir2;
  set pb_per_ir;
  if count > 1;
run; *good: 1 match per ireckon;

proc sort data=ir2pacbio;
   by sample_id ir_transcript_id;
proc sort data=ireckon_refseq_best;
   by sample_id ir_transcript_id;
run;

data ireckon_pacbio_best;
  merge ireckon_refseq_best (in=in1) ir2pacbio (in=in2);
  by sample_id ir_transcript_id;
  if in2 then flag_has_pacbio_hit=1; else flag_has_pacbio_hit=0;
run;

/* Reduced Ref 1 transcripts */

data xs100;
  set event.bin_xscripts_by_dtct_apn0;
  where perc_features_dtct=1;
  keep transcript_id;
run;

proc sort data=ireckon_pacbio_best;
  by transcript_id;
proc sort data=xs100;
  by transcript_id;
run;

data ireckon_rs_pb_rr1;
  merge ireckon_pacbio_best (in=in1) xs100 (in=in2);
  by transcript_id;
  if in2 then flag_in_xs100_apn0=1; else flag_in_xs100_apn0=0;
  if in1 then output;
run;

/* Reduced Ref 2 transcripts */

data xs75;
  set event.bin_xscripts_by_dtct_apn5;
  where perc_features_dtct ge 0.75;
  keep transcript_id;
run;

proc sort data=ireckon_rs_pb_rr1;
  by transcript_id;
proc sort data=xs75;
  by transcript_id;
run;

data ireckon_rs_pb_rr2;
  merge ireckon_rs_pb_rr1 (in=in1) xs75 (in=in2);
  by transcript_id;
  if in2 then flag_in_xs75_apn5=1; else flag_in_xs75_apn5=0;
  if in1 then output;
run;

/* Now I want to split by sample, and figure out intersections, unions, etc.
I want to also add in the "transcript type" if possible
*/

data ir_npc1 ir_npc2;
  set ireckon_rs_pb_rr2;
  if sample_id="NSC1" then output ir_npc1;
  else if sample_id="NSC2" then output ir_npc2;
run;

data ir_npc1_2;
   set ir_npc1;
   keep ir_transcript_id flag_has_refseq_hit flag_has_pacbio_hit
        flag_in_xs100_apn0 flag_in_xs75_apn5;
   rename flag_has_refseq_hit=flag_npc1_refseq_hit
          flag_has_pacbio_hit=flag_npc1_pacbio_hit
          flag_in_xs100_apn0=flag_npc1_in_xs100
          flag_in_xs75_apn5=flag_npc1_in_xs75;
run;

data ir_npc2_2;
   set ir_npc2;
   keep ir_transcript_id flag_has_refseq_hit flag_has_pacbio_hit
        flag_in_xs100_apn0 flag_in_xs75_apn5;
   rename flag_has_refseq_hit=flag_npc2_refseq_hit
          flag_has_pacbio_hit=flag_npc2_pacbio_hit
          flag_in_xs100_apn0=flag_npc2_in_xs100
          flag_in_xs75_apn5=flag_npc2_in_xs75;
run;

proc sort data=ir_npc1_2;
  by ir_transcript_id;
proc sort data=ir_npc2_2;
  by ir_transcript_id;
run;

proc means data=ir_npc1_2 noprint;
   by ir_transcript_id;
   var flag_npc1_refseq_hit flag_npc1_pacbio_hit flag_npc1_in_xs100 flag_npc1_in_xs75;
   output out=ir_npc1_3(drop=_TYPE_ _FREQ_) max=;
run;

proc means data=ir_npc2_2 noprint;
   by ir_transcript_id;
   var flag_npc2_refseq_hit flag_npc2_pacbio_hit flag_npc2_in_xs100 flag_npc2_in_xs75;
   output out=ir_npc2_3(drop=_TYPE_ _FREQ_) max=;
run;


data ir_npc1_npc2;
   merge ir_npc1_3 (in=in1) ir_npc2_3 (in=in2);
   by ir_transcript_id;
   if in1 then flag_in_npc1=1; else flag_in_npc1=0;
   if in2 then flag_in_npc2=1; else flag_in_npc2=0;
run;

data ir_npc1_npc2_2;
  set ir_npc1_npc2;
  array change _numeric_;
  do over change;
    if change=. then change=0;
    end;
run;

/* Add in flags for looking at the intersection and union 
   If there is a PB or Refseq hit in one of the reps, then there is a hit for the
   union. For the intersection, transcript must be detected in both */


data ir_npc1_npc2_3;
  set ir_npc1_npc2_2;
  /* Union */
  if flag_in_npc1=1 or flag_in_npc2=1 then flag_in_npc_any=1;
  else flag_in_npc_any=1;
  if flag_npc1_refseq_hit=1 or flag_npc2_refseq_hit=1 then flag_any_refseq_hit=1;
  else flag_any_refseq_hit=0;
  if flag_npc1_pacbio_hit=1 or flag_npc2_pacbio_hit=1 then flag_any_pacbio_hit=1;
  else flag_any_pacbio_hit=0;
  if flag_npc1_in_xs100=1 or flag_npc2_in_xs100=1 then flag_any_in_xs100=1;
  else flag_any_in_xs100=0;
  if flag_npc1_in_xs75=1 or flag_npc2_in_xs75=1 then flag_any_in_xs75=1;
  else flag_any_in_xs75=0;
  /* Intersection -- here the transcript has to be in both reps */
  if flag_in_npc1=1 and flag_in_npc2=1 then do;
      flag_in_npc_all=1;
      if flag_npc1_refseq_hit=1 or flag_npc2_refseq_hit=1 then flag_all_refseq_hit=1;
      else flag_all_refseq_hit=0;
      if flag_npc1_pacbio_hit=1 or flag_npc2_pacbio_hit=1 then flag_all_pacbio_hit=1;
      else flag_all_pacbio_hit=0;
      if flag_npc1_in_xs100=1 or flag_npc2_in_xs100=1 then flag_all_in_xs100=1;
      else flag_all_in_xs100=0;
      if flag_npc1_in_xs75=1 or flag_npc2_in_xs75=1 then flag_all_in_xs75=1;
      else flag_all_in_xs75=0;
      end;
  else do;
      flag_in_npc_all=0;
      flag_all_refseq_hit=0;
      flag_all_pacbio_hit=0;
      flag_all_in_xs100=0;
      flag_all_in_xs75=0;
      end;
run;


/* Count overlap between iReckon, RefSeq, Pacbio, and reduced references */

proc freq data=ir_npc1_npc2_3 noprint;
   /* NPC1 */
   tables flag_in_npc1*flag_npc1_refseq_hit*flag_npc1_pacbio_hit*
          flag_npc1_in_xs100*flag_npc1_in_xs75 / out=ir_npc1_counts;
   /* NPC2 */
   tables flag_in_npc2*flag_npc2_refseq_hit*flag_npc2_pacbio_hit*
          flag_npc2_in_xs100*flag_npc2_in_xs75 / out=ir_npc2_counts;
   /* Union */
   tables flag_in_npc_any*flag_any_refseq_hit*flag_any_pacbio_hit*
          flag_any_in_xs100*flag_any_in_xs75 / out=ir_union_counts;
   /* Intersection */
   tables flag_in_npc_all*flag_all_refseq_hit*flag_all_pacbio_hit*
          flag_all_in_xs100*flag_all_in_xs75 / out=ir_intersect_counts;
run;

proc print data=ir_npc1_counts;
run;
proc print data=ir_npc2_counts;
run;
proc print data=ir_union_counts;
run;
proc print data=ir_intersect_counts;
run;


/*

                                          flag_       flag_
flag_in_    flag_npc1_    flag_npc1_    npc1_in_    npc1_in_
  npc1      refseq_hit    pacbio_hit      xs100       xs75      COUNT

    0            0             0            0           0        9120
    1            0             0            0           0       11630
    1            0             1            0           0        1278
    1            1             0            0           0         355
    1            1             0            0           1          50
    1            1             0            1           0         124
    1            1             0            1           1         172
    1            1             1            0           0         256
    1            1             1            0           1          62
    1            1             1            1           0          36
    1            1             1            1           1         355

                                           flag_       flag_
 flag_in_    flag_npc2_    flag_npc2_    npc2_in_    npc2_in_
   npc2      refseq_hit    pacbio_hit      xs100       xs75      COUNT

     0            0             0            0           0       11856
     1            0             0            0           0        9675
     1            0             1            0           0         825
     1            1             0            0           0         264
     1            1             0            0           1          52
     1            1             0            1           0          68
     1            1             0            1           1         167
     1            1             1            0           0         181
     1            1             1            0           1          70
     1            1             1            1           0          14
     1            1             1            1           1         266


                                           flag_      flag_
  flag_in_     flag_any_     flag_any_    any_in_    any_in_
   npc_any    refseq_hit    pacbio_hit     xs100       xs75     COUNT

      1            0             0           0          0       19164
      1            0             1           0          0        1948
      1            1             0           0          0         586
      1            1             0           0          1          97
      1            1             0           1          0         160
      1            1             0           1          1         323
      1            1             1           0          0         411
      1            1             1           0          1         126
      1            1             1           1          0          48
      1            1             1           1          1         575


                                          flag_      flag_
 flag_in_     flag_all_     flag_all_    all_in_    all_in_
  npc_all    refseq_hit    pacbio_hit     xs100       xs75     COUNT

     0            0             0           0          0       20976
     1            0             0           0          0        2108
     1            0             1           0          0         160
     1            1             0           0          0          38
     1            1             0           0          1           6
     1            1             0           1          0          32
     1            1             0           1          1          22
     1            1             1           0          0          34
     1            1             1           0          1           6
     1            1             1           1          0           2
     1            1             1           1          1          54

*/


