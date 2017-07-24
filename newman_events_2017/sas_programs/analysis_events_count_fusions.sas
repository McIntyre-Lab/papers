ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of fusions (exon regions) that are expressed at APN>0 */

data fusions;
  set event.flag_fusion_on;
  keep fusion_id flag_fusion_nsc_on;
  rename flag_fusion_nsc_on=flag_fusion_on;
run;


proc freq data=fusions;
   tables flag_fusion_on;
run;


/*
    flag_fusion_                             Cumulative    Cumulative
          on    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       96413       35.04         96413        35.04
               1      178764       64.96        275177       100.00


*/


data event.fusions_on_apn_gt0;
   set fusions;
run;

