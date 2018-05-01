ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* If an unannotated junction has at least one hit, then flag the junction has having Pacbio support
   and look at the cross-tabulation with IR classification  */

data blast_hits;
 set event.unannot_junc_best_blast_hits;
 if perc_identity < 0.95 then delete;
run;

data blast_hits_uniq;
  set blast_hits;
  keep event_id;
run;

proc sort data=ir_hits_uniq nodup;
   by event_id;
run;

/* Merge with full set of classified IR events and look at cross-over */

data junc_blasted;
   set event.unannot_junc_pacbio_hits;
   keep event_id;
run;


proc sort data=junc_blasted nodup;
  by event_id;
proc sort data=blast_hits_uniq nodup;
  by event_id;
run;

data junc_flag_hits;
   merge junc_blasted (in=in1) blast_hits_uniq (in=in2);
   by event_id;
   if in1 then flag_ir_blast_hit=1; else flag_ir_blast_hit=0;
   if in2 then flag_ir_good_hit=1; else flag_ir_good_hit=0;
run;

proc freq data=junc_flag_hits noprint;
   tables flag_ir_blast_hit*flag_ir_good_hit / out=junc_counts;
proc print data=junc_counts;
run;


/*
 
SUMMARY:
  flag_ir_    flag_ir_
 blast_hit    good_hit    COUNT

     1            0         164
     1            1        1399


1399 junctions with a good blast hit
164 without

*/

/* Now I want to count by PB transcript type
  
   */

proc sort data=blast_hits;
  by event_id;
data blast_hits2class;
  set blast_hits;
  keep event_id pb_status;
run;

proc sort data=blast_hits2class nodup;
   by event_id pb_status;
proc freq data=blast_hits2class noprint;
   tables event_id / out=pb_count;
proc sort data=pb_count;
  by descending count;
run; *max 4;


data cat_pb; 
  array pb[2] $ 30;
  retain pb1-pb2;
  set blast_hits2class;
  by event_id;
  if first.event_id then do;
     call missing(of pb1-pb2);
     records = 0;
  end;
  records + 1;
  pb[records]=pb_status;
  if last.event_id then output;
run;

  *clean up the output file;
data cat_pb2;
  set cat_pb;
  length pb_hit_type $ 150;
         pb_hit_type= catx("|", OF pb1-pb2);
  keep event_id pb_hit_type;
  run;

/* Merge in IR class flags */

proc freq data=cat_pb2;
   tables pb_hit_type;
run;


/*
 pb_hit_type    Frequency     Percent
 --------------------------------------
 Known               796       56.90
 Known|Novel         260       18.58
 Novel               343       24.52

*/

