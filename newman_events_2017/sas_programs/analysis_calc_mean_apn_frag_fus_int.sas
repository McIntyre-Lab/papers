ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For fusions, fragments, introns, calc the mean APN (NPCs for now) */

data frag_counts;
  set event.mm10_refseq_fragment_counts;
  where sample_id in ('NSC1','NSC2');
  keep sample_id fragment_id apn;
run;

data fus_counts;
  set event.mm10_refseq_fusion_counts;
  where sample_id in ('NSC1','NSC2');
  keep sample_id fusion_id apn;
run;

data int_counts;
  set event.mm10_refseq_intron_counts;
  where sample_id in ('NSC1','NSC2');
  keep sample_id intron_id apn;
run;


data ir_counts;
   set refseq.rfsq_counts_by_splicing;
   where event_id ? "intron" and sample_id in ('NSC1','NSC2');
run;

proc sort data=frag_counts;
   by fragment_id;
proc means data=frag_counts noprint;
   by fragment_id;
   var apn;
   output out=mean_frag_apn mean=;
run;

proc sort data=fus_counts;
   by fusion_id;
proc means data=fus_counts noprint;
   by fusion_id;
   var apn;
   output out=mean_fus_apn mean=;
run;

proc sort data=int_counts;
   by intron_id;
proc means data=int_counts noprint;
   by intron_id;
   var apn;
   output out=mean_int_apn mean=;
run;


proc sort data=ir_counts;
   by event_id;
proc means data=ir_counts noprint;
   by event_id;
   var apn;
   output out=mean_ir_apn mean=;
run;

/* Make permenant */

data event.mean_apn_intron_nsc;
  set mean_int_apn;
  keep intron_id apn;
  rename apn=mean_apn_intron;
run;

data event.mean_apn_fragment_nsc;
  set mean_frag_apn;
  keep fragment_id apn;
  rename apn=mean_apn_fragment;
run;

data event.mean_apn_fusion_nsc;
  set mean_fus_apn;
  keep fusion_id apn;
  rename apn=mean_apn_fusion;
run;

data event.mean_apn_IR_nsc;
  set mean_ir_apn;
  keep event_id apn;
  rename apn=mean_apn_IR;
run;


