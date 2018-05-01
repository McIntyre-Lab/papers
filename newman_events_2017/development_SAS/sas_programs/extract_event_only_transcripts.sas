/* need to make PacBio-plus transcriptomes, so first I need to take the Refseq lists and exclude ones with PB matches */

ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';

/* get APN<0/5 50% and 75% lists */

data list1 list2;
  set event.xscripts_w_unique_by_bin;
  if perc_features_dtct >= 0.5 then output list1;
  if perc_features_dtct >= 0.75 then output list2;
   keep transcript_id;
run;

data list3 list4;
  set event.bin_xscripts_by_dtct_apn5;
  if perc_features_dtct >= 0.5 then output list3;
  if perc_features_dtct >= 0.75 then output list4;
   keep transcript_id;
run;

 
  
data pb_gene2xs;
   set event.pacbio_transcripts;
   keep pacbio_gene_id pacbio_id;
run;

     data WORK.SPLICE_MATCH    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/references/pacbio_isoforms_list_for_import.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat pacbio_id $10. ; informat source $6.;    informat feature_type $10. ;
        informat start best32. ;  informat end best32.;   informat frame best32. ;
        informat strand $1. ;     informat score best32.; informat transcript_id $18. ;
        informat match_type $23.; informat note $28. ;
        format pacbio_id $10. ;   format source $6. ;     format feature_type $10. ;
        format start best12. ;    format end best12. ;    format frame best12. ;
        format strand $1. ;       format score best12. ;  format transcript_id $18. ;
        format match_type $23. ;  format note $28. ;
        input pacbio_id $ source $ feature_type $ start end
              frame strand $ score transcript_id $ match_type $ note $
              ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data pb2refseq;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep pacbio_id splice_match_id ;
   rename splice_match_id=transcript_id;
run; *8722;


data pb2keep;
  set event.pacbio2refseq_id_nomulti;
  keep pacbio_id;
run;

proc sort data=pb2keep nodup;
  by pacbio_id;
proc sort data=pb2refseq nodup;
  by pacbio_id transcript_id;
run;

data xs_w_pb;
  merge pb2refseq (in=in1) pb2keep (in=in2);
  by pacbio_id;
  if in1 and in2;
run;


proc sort data=xs_w_pb;
   by transcript_id;
proc sort data=list1;
   by transcript_id;
proc sort data=xs_w_pb;
   by transcript_id;
proc sort data=list2;
   by transcript_id;
proc sort data=xs_w_pb;
   by transcript_id;
proc sort data=list3;
   by transcript_id;
proc sort data=xs_w_pb;
   by transcript_id;
proc sort data=list4;
   by transcript_id;
run;



data list1_2;
   merge list1 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   if in1 and in2 then delete;
   if in1 then output;
   keep transcript_id;
   run;
   
   data list2_2;
   merge list2 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   if in1 and in2 then delete;
   if in1 then output;
      keep transcript_id;
   run;
   
   data list3_2;
   merge list3 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   if in1 and in2 then delete;
   if in1 then output;
      keep transcript_id;
   run;
   
   data list4_2;
   merge list4 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   if in1 and in2 then delete;
   if in1 then output;
      keep transcript_id;
   run;
   
proc export data=list1_2
outfile="!MCLAB/event_analysis/references/event_only_apn0_50perc.txt"
dbms=tab replace; putnames=no;
run;

proc export data=list2_2
outfile="!MCLAB/event_analysis/references/event_only_apn0_75perc.txt"
dbms=tab replace; putnames=no;
run;

proc export data=list3_2
outfile="!MCLAB/event_analysis/references/event_only_apn5_50perc.txt"
dbms=tab replace; putnames=no;
run;

proc export data=list4_2
outfile="!MCLAB/event_analysis/references/event_only_apn5_75perc.txt"
dbms=tab replace; putnames=no;
run;

