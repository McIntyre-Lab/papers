
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
   keep transcript_id pacbio_id;
   run;

   data list2_2;
   merge list2 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   keep transcript_id pacbio_id;
   run;

   data list3_2;
   merge list3 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
  keep transcript_id pacbio_id;
   run;

   data list4_2;
  merge list4 (in=in1) xs_w_pb (in=in2);
   by transcript_id;
   keep transcript_id pacbio_id ;
   run;

data pb_all;
   set event.pacbio_transcripts;
   keep pacbio_gene_id pacbio_id;
run;

proc sort data=list1_2;
  by pacbio_id;
  proc sort data=list2_2;
  by pacbio_id;
  proc sort data=list3_2;
  by pacbio_id;
  proc sort data=list4_2;
  by pacbio_id;
  proc sort data=pb_all;
  by pacbio_id;
  run;

data list1_3;
   merge list1_2 pb_all;
   by pacbio_id;
run;
data list2_3;
   merge list2_2 pb_all;
   by pacbio_id;
run;
data list3_3;
   merge list3_2 pb_all;
   by pacbio_id;
run;
data list4_3;
   merge list4_2 pb_all;
   by pacbio_id;
run;



data xs2gene;
  set event.feature2xs2gene_exp_only_nomulti;
  keep transcript_id gene_id;
  run;

proc sort data=xs2gene nodup;
  by transcript_id gene_id;
proc sort data=list1_3;
  by transcript_id;
proc sort data=list2_3;
  by transcript_id;
proc sort data=list3_3;
  by transcript_id;
proc sort data=list4_3;
  by transcript_id;
  run;
  
  data list1_w_gene;
     merge xs2gene (in=in1) list1_3 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;
     
  data list2_w_gene;
     merge xs2gene (in=in1) list2_3 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;
     
  data list3_w_gene;
     merge xs2gene (in=in1) list3_3 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;
     
  data list4_w_gene;
     merge xs2gene (in=in1) list4_3 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;
  
 
data pb2gene;
   set event.pacbio2refseq_gene_nomulti;
   keep gene_id pacbio_gene_id;
   rename gene_id=refseq_gene_id;
   run;
   

 proc sort data=list1_w_gene;
    by pacbio_gene_id;
  proc sort data=list2_w_gene;
    by pacbio_gene_id;
  proc sort data=list3_w_gene;
    by pacbio_gene_id;
 proc sort data=list4_w_gene;
   by pacbio_gene_id;
 proc sort data=pb2gene nodup;
   by pacbio_gene_id;
  run;
  
  data pb2gene_1;
     merge list1_w_gene (in=in1) pb2gene (in=in2);
     by pacbio_gene_id;
     if in1 then output;
  run;
    data pb2gene_2;
     merge list2_w_gene (in=in1) pb2gene (in=in2);
     by pacbio_gene_id;
     if in1 then output;
  run;
    data pb2gene_3;
     merge list3_w_gene (in=in1) pb2gene (in=in2);
     by pacbio_gene_id;
     if in1 then output;
  run;
    data pb2gene_4;
     merge list4_w_gene (in=in1) pb2gene (in=in2);
     by pacbio_gene_id;
     if in1 then output;
  run;
  
  data pb2gene_1a;
     length gene_id2 $15.;
     length transcript_id2 $15.;
     set pb2gene_1;
     if gene_id ne '' then gene_id2=gene_id;
     else if pacbio_gene_id ne '' then gene_id2=pacbio_gene_id;
     else gene_id2="NovelGene";
     
     if pacbio_id ne '' then transcript_id2=pacbio_id;
     else transcript_id2=transcript_id;
     keep gene_id2 transcript_id2;
     run;
     
       data pb2gene_2a;
     length gene_id2 $15.;
     length transcript_id2 $15.;
     set pb2gene_2;
     if gene_id ne '' then gene_id2=gene_id;
     else if pacbio_gene_id ne '' then gene_id2=pacbio_gene_id;
     else gene_id2="NovelGene";
     
     if pacbio_id ne '' then transcript_id2=pacbio_id;
     else transcript_id2=transcript_id;
     keep gene_id2 transcript_id2;
     run;
     
       data pb2gene_3a;
     length gene_id2 $15.;
     length transcript_id2 $15.;
     set pb2gene_3;
     if gene_id ne '' then gene_id2=gene_id;
     else if pacbio_gene_id ne '' then gene_id2=pacbio_gene_id;
     else gene_id2="NovelGene";
     
     if pacbio_id ne '' then transcript_id2=pacbio_id;
     else transcript_id2=transcript_id;
     keep gene_id2 transcript_id2;
     run;
     
       data pb2gene_4a;
     length gene_id2 $15.;
     length transcript_id2 $15.;
     set pb2gene_4;
     if gene_id ne '' then gene_id2=gene_id;
     else if pacbio_gene_id ne '' then gene_id2=pacbio_gene_id;
     else gene_id2="NovelGene";
     
     if pacbio_id ne '' then transcript_id2=pacbio_id;
     else transcript_id2=transcript_id;
     keep gene_id2 transcript_id2;
     run;
     
     

  
    proc export data=pb2gene_1a
       outfile="!MCLAB/event_analysis/pacbio_plus_apn0_50perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
       
    proc export data=pb2gene_2a
       outfile="!MCLAB/event_analysis/pacbio_plus_apn0_75perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
       
    proc export data=pb2gene_3a
       outfile="!MCLAB/event_analysis/pacbio_plus_apn5_50perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
       
    proc export data=pb2gene_4a
       outfile="!MCLAB/event_analysis/pacbio_plus_apn5_75perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
       
