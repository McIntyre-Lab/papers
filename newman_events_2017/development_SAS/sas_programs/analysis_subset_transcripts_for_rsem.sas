libname event '!MCLAB/event_analysis/sas_data';

/* Export list of transcripts for RSEM 

I want to make the following lists:

1. All Refseq (so just use original input)
2. All Pacbio (again, raw input)
3. all 73K transcripts detected
4. exp transcripts with all features detected at APN>0
5. exp transcripts with at least 75% features detected at APN>0

*/


data xscript1 xscript2 xscript3;
  set event.xscripts_w_unique_by_bin;
  if perc_features_dtct > 0 then output xscript1;
  if perc_features_dtct >= 0.75 then output xscript2;
  if perc_features_dtct = 1 then output xscript3;
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

/* 75% detected, APN 5 and 10 */

data xscript4;
   set event.bin_xscripts_by_dtct_apn5;
   where perc_features_dtct ge 0.75;
   keep transcript_id;
run;

data xscript5;
   set event.bin_xscripts_by_dtct_apn10;
   where perc_features_dtct ge 0.75;
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

/* Now for each transcript list, I need to make a transcript-to-gene index where the first column is geneID
   and the second is transcriptID. This is used by RSEM to calculate gene-level expression. While we don't need
   this, I want to include it in case RSEM does something funky with transcripts of the same gene */

data rs_gene2xs;
   retain gene_id transcript_id;
   set event.feature2xs2gene;
   keep gene_id transcript_id;
run;

data pb_gene2xs;
   set event.pacbio_transcripts;
   keep pacbio_gene_id pacbio_id;
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


proc sort data=pb_gene2xs nodup;
   by pacbio_id pacbio_gene_id;
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



proc export data=pb_gene2xs outfile="!MCLAB/event_analysis/references/pacbio_mm10_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;






