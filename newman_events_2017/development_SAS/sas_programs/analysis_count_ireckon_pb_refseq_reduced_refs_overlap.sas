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
   if log_rpkm=0 then rpkm_bin=0;
   else if log_rpkm < 0.1 then rpkm_bin=1;
   else if log_rpkm < 0.5 then rpkm_bin=2;
   else if log_rpkm < 1.5 then rpkm_bin=3;
   else rpkm_bin=4;
   keep ir_transcript_id sample_id rpkm_bin;
run;

proc sort data=ireckon;
   by sample_id ir_transcript_id;
proc means data=ireckon noprint;
   by sample_id ir_transcript_id;
   var rpkm_bin;
   output out=ireckon2(drop=_TYPE_ _FREQ_) max=;
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
proc sort data=ireckon2;
   by sample_id ir_transcript_id;
run;

data ireckon2_refseq_best;
  merge ireckon2 (in=in1) ir2refseq (in=in2);
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
proc sort data=ireckon2_refseq_best;
   by sample_id ir_transcript_id;
run;

data ireckon2_pacbio_best;
  merge ireckon2_refseq_best (in=in1) ir2pacbio (in=in2);
  by sample_id ir_transcript_id;
  if in2 then flag_has_pacbio_hit=1; else flag_has_pacbio_hit=0;
run;

/* Reduced Ref 1 transcripts */

data xs100;
  set event.bin_xscripts_by_dtct_apn0;
  where perc_features_dtct=1;
  keep transcript_id;
run;

proc sort data=ireckon2_pacbio_best;
  by transcript_id;
proc sort data=xs100;
  by transcript_id;
run;

data ireckon_rs_pb_rr1;
  merge ireckon2_pacbio_best (in=in1) xs100 (in=in2);
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

/* Count overlap between iReckon, RefSeq, Pacbio, and reduced references for each sample */

proc sort data=ireckon_rs_pb_rr2;
  by sample_id rpkm_bin;
proc freq data=ireckon_rs_pb_rr2 noprint;
  by sample_id;
  tables flag_has_refseq_hit*flag_has_pacbio_hit*flag_in_xs100_apn0*flag_in_xs75_apn5 / out=xs_count_by_sample;
run;

proc print data=xs_count_by_sample;
  where sample_id="NSC1";
run;

/*
                                          flag_in_    flag_in_
   sample_     flag_has_     flag_has_     xs100_       xs75_
     id       refseq_hit    pacbio_hit      apn0        apn5      COUNT

    NSC1           0             0            0           0       11804
    NSC1           0             1            0           0        1104
    NSC1           1             0            0           0         364
    NSC1           1             0            0           1          56
    NSC1           1             0            1           0         150
    NSC1           1             0            1           1         203
    NSC1           1             1            0           0         247
    NSC1           1             1            0           1          56
    NSC1           1             1            1           0          10
    NSC1           1             1            1           1         324

*/

/* Export */

proc export data=xs_count_by_sample
   outfile="!MCLAB/event_analysis/analysis_output/iReckon_vs_Refseq_PacBio_Events_by_sample.csv"
   dbms=csv replace;
run;

