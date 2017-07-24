ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* Export list of transcripts for RSEM 

I want to make the following lists:

1. All Refseq (so just use original input)
2. All Pacbio (again, raw input)
3. all 73K transcripts detected
4. exp transcripts with all features detected at APN>0
5. exp transcripts with at least 75% features detected at APN>0

*/


data xscript1 xscript2 xscript3 xscript50;
  set event.xscripts_w_unique_by_bin;
  if perc_features_dtct > 0 then output xscript1;
  if perc_features_dtct >= 0.75 then output xscript2;
  if perc_features_dtct = 1 then output xscript3;
  if perc_features_dtct >= 0.5 then output xscript50;
  keep transcript_id;
run;

proc export data=xscript1 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts.txt"
   dbms=tab replace;
   putnames=no;
run;

proc export data=xscript2 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_75perc_dtct.txt"
   dbms=tab replace;
   putnames=no;
run;

proc export data=xscript3 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_100perc_dtct.txt"
   dbms=tab replace;
   putnames=no;
run;


proc export data=xscript50 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_50perc_dtct.txt"
   dbms=tab replace;
   putnames=no;
run;

/* 75% detected, APN 5 and 10 */

data xscript4 xscript4_100 xscript4_50 xscript_apn5_for_pb;
   set event.bin_xscripts_by_dtct_apn5;
  if perc_features_dtct > 0 then output xscript_apn5_for_pb;
  if perc_features_dtct >= 0.75 then output xscript4;
  if perc_features_dtct = 1 then output xscript4_100;
  if perc_features_dtct >= 0.5 then output xscript4_50;
   keep transcript_id;
run;

data xscript5 xscript5_100 xscript5_50 xscript_apn10_for_pb;
   set event.bin_xscripts_by_dtct_apn10;
  if perc_features_dtct > 0 then output xscript_apn10_for_pb;
  if perc_features_dtct >= 0.75 then output xscript5;
  if perc_features_dtct = 1 then output xscript5_100;
  if perc_features_dtct >= 0.5 then output xscript5_50;
   keep transcript_id;
run;


proc export data=xscript4 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_75perc_dtct_apn5.txt"
   dbms=tab replace;
   putnames=no;
run;

proc export data=xscript5 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_75perc_dtct_apn10.txt"
   dbms=tab replace;
   putnames=no;
run;


proc export data=xscript4_100 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_100perc_dtct_apn5.txt"
   dbms=tab replace;
   putnames=no;
run;

proc export data=xscript5_100 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_100perc_dtct_apn10.txt"
   dbms=tab replace;
   putnames=no;
run;


proc export data=xscript4_50 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_50perc_dtct_apn5.txt"
   dbms=tab replace;
   putnames=no;
run;

proc export data=xscript5_50 outfile="!MCLAB/event_analysis/analysis_output/refseq_list_exp_transcripts_50perc_dtct_apn10.txt"
   dbms=tab replace;
   putnames=no;
run;


/* Now for each transcript list, I need to make a transcript-to-gene index where the first column is geneID
   and the second is transcriptID. This is used by RSEM to calculate gene-level expression. While we don't need
   this, I want to include it in case RSEM does something funky with transcripts of the same gene */

data rs_gene2xs;
   retain gene_id transcript_id;
   set event.feature2xs2gene;
   keep gene_id transcript_id;
run;

proc sort data=rs_gene2xs nodup;
   by transcript_id gene_id;
proc sort data=xscript1 nodup;
   by transcript_id;
proc sort data=xscript2 nodup;
   by transcript_id ;
proc sort data=xscript3 nodup;
   by transcript_id ;
proc sort data=xscript4 nodup;
   by transcript_id ;
proc sort data=xscript5 nodup;
   by transcript_id ;
proc sort data=xscript4_100 nodup;
   by transcript_id ;
proc sort data=xscript5_100 nodup;
   by transcript_id ;
proc sort data=xscript4_50 nodup;
   by transcript_id ;
proc sort data=xscript5_50 nodup;
   by transcript_id ;
proc sort data=xscript50 nodup;
   by transcript_id ;
run;

data gene2xs_list1;
  merge rs_gene2xs (in=in1) xscript1 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list2;
  merge rs_gene2xs (in=in1) xscript2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list3;
  merge rs_gene2xs (in=in1) xscript3 (in=in2);
  by transcript_id;
  if in1 and in2;
run;


data gene2xs_list4;
  merge rs_gene2xs (in=in1) xscript4 (in=in2);
  by transcript_id;
  if in1 and in2;
run;



data gene2xs_list5;
  merge rs_gene2xs (in=in1) xscript5 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list4_100;
  merge rs_gene2xs (in=in1) xscript4_100 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list5_100;
  merge rs_gene2xs (in=in1) xscript5_100 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list4_50;
  merge rs_gene2xs (in=in1) xscript4_50 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list5_50;
  merge rs_gene2xs (in=in1) xscript5_50 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data gene2xs_list50;
  merge rs_gene2xs (in=in1) xscript50 (in=in2);
  by transcript_id;
  if in1 and in2;
run;


/* Export lists */

proc export data=rs_gene2xs outfile="!MCLAB/event_analysis/references/refseq_mm10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list1 outfile="!MCLAB/event_analysis/references/refseq_list_exp_transcripts_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list2 outfile="!MCLAB/event_analysis/references/refseq_list_exp_transcripts_75perc_dtct_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list3 outfile="!MCLAB/event_analysis/references/refseq_list_exp_transcripts_100perc_dtct_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;


proc export data=gene2xs_list4 outfile="!MCLAB/event_analysis/references/refseq_list_exp_transcripts_75perc_dtct_apn5_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;


proc export data=gene2xs_list5 outfile="!MCLAB/event_analysis/references/refseq_list_exp_transcripts_75perc_dtct_apn10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;





proc export data=gene2xs_list50 outfile="!MCLAB/event_analysis/references/refseq_mm10_exp_transcripts_50perc_dtct_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list4_100 outfile="!MCLAB/event_analysis/references/refseq_mm10_exp_transcripts_100perc_dtct_apn5_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list5_100 outfile="!MCLAB/event_analysis/references/refseq_mm10_exp_transcripts_100perc_dtct_apn10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list4_50 outfile="!MCLAB/event_analysis/references/refseq_mm10_exp_transcripts_50perc_dtct_apn5_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list5_50 outfile="!MCLAB/event_analysis/references/refseq_mm10_exp_transcripts_50perc_dtct_apn10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;




/* PACBIO */

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
proc sort data=xscript1;
   by transcript_id;
proc sort data=xscript_apn5_for_pb;
   by transcript_id;
proc sort data=xscript_apn10_for_pb;
   by transcript_id;
run;

data pb_list1;
  merge xs_w_pb (in=in1) xscript1 (in=in2);
  by transcript_id;
  if in1 and in2;
  keep pacbio_id;
run;


data pb_list2;
  merge xs_w_pb (in=in1) xscript_apn5_for_pb (in=in2);
  by transcript_id;
  if in1 and in2;
  keep pacbio_id;
run;


data pb_list3;
  merge xs_w_pb (in=in1) xscript_apn10_for_pb (in=in2);
  by transcript_id;
  if in1 and in2;
  keep pacbio_id;
run;

proc sort data=pb_gene2xs nodup;
   by pacbio_id pacbio_gene_id;
proc sort data=pb_list1 nodup;
   by pacbio_id;
proc sort data=pb_list2 nodup;
   by pacbio_id;
proc sort data=pb_list3 nodup;
   by pacbio_id ;
run;

data pb_gene2xs_list1;
  merge pb_gene2xs (in=in1) pb_list1 (in=in2);
  by pacbio_id;
  if in1 and in2;
run;

data pb_gene2xs_list2;
  merge pb_gene2xs (in=in1) pb_list2 (in=in2);
  by pacbio_id;
  if in1 and in2;
run;

data pb_gene2xs_list3;
  merge pb_gene2xs (in=in1) pb_list3 (in=in2);
  by pacbio_id;
  if in1 and in2;
run;



proc export data=pb_gene2xs outfile="!MCLAB/event_analysis/references/pacbio_mm10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=pb_gene2xs_list1 outfile="!MCLAB/event_analysis/references/pacbio_mm10_apn0_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=pb_gene2xs_list2 outfile="!MCLAB/event_analysis/references/pacbio_mm10_apn5_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=pb_gene2xs_list3 outfile="!MCLAB/event_analysis/references/pacbio_mm10_apn10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;




proc export data=pb_list1 outfile="!MCLAB/event_analysis/analysis_output/pacbio_mm10_apn0.txt"
   dbms=tab replace;
   putnames=no;
run;
proc export data=pb_list2 outfile="!MCLAB/event_analysis/analysis_output/pacbio_mm10_apn5.txt"
   dbms=tab replace;
   putnames=no;
run;
proc export data=pb_list3 outfile="!MCLAB/event_analysis/analysis_output/pacbio_mm10_apn10.txt"
   dbms=tab replace;
   putnames=no;
run;





