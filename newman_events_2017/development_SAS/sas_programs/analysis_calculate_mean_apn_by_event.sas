ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Calculate mean APN for junctions and fragments */

data frag_counts;
    set event.mm10_refseq_fragment_counts;
    where sample_id ? "NSC";
    keep sample_id fragment_id apn;
    rename fragment_id=feature_id;
run;

data event_counts;
    set event.mm10_refseq_splicing_counts;
    where sample_id ? "NSC";
    keep sample_id event_id apn;
    rename event_id=feature_id;
run;

data feature_counts;
   set event_counts frag_counts;
run;

proc sort data=feature_counts;
   by feature_id sample_id;
proc means data=feature_counts noprint;
   by feature_id;
   var apn;
   output out=mean_apn_by_feat mean=mean_apn_npc;
run;

/* Make permenant */

data event.mean_apn_events_npc;
   set mean_apn_by_feat;
   drop _TYPE_ _FREQ_;
run;

