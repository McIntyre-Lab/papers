ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


data hits_annot;
   set event.blast_junc_to_pb_w_annot;
   where flag_junction_annotated=0 and flag_feature_on_ge5=1;
run;

data ir2keep;
  set event.ir_reclassification_v3_apn5;
  where flag_splicing_on_ge5=1 and flag_low_expressed=0;
  keep event_id;
run;

data junc2keep;
   set event.unannot_junc_dtct_flag_fus;
   keep event_id;
run;

data mean_apn;
   set event.mean_apn_events_npc;
   rename feature_id=event_id;
run;

proc sort data=hits_annot;
   by event_id;
proc sort data=mean_apn;
   by event_id;
proc sort data=ir2keep;
   by event_id;
proc sort data=junc2keep;
   by event_id;
run;

data hits_w_apn ;
  merge hits_annot (in=in1) mean_apn (in=in2) ir2keep (in=in3) junc2keep (in=in4);
  by event_id;
  if in1 and in2 and in3 then output hits_w_apn;
  if in1 and in2 and in4 then output hits_w_apn;
run;


proc export data=hits_w_apn
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_hits_unannotated_events_apn.csv"
     dbms=csv replace;
run;

