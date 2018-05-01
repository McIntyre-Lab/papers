/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Of the analyzable PacBio transcripts, what is the detection rate in terms of novel junctions/isoforms
   for STAR, Event Analysis, and iReckon?

   Transcript-based:
   (1) Identify the set of novel PacBio transcripts
   (2) What is the overlap between iReckon transcripts and novel PacBio?
		-> NSC1, NSC2, NSC1 ∪ NSC2, NSC1 ∩ NSC2
		-> Cross-tabulation
	(note this will be more-or-less the same for the general iReckon vs PacBio comparison
     so code here can be used again)
*/

/* Order of operations:
(1) iReckon transcript BLAST table: NSC1, NSC2, NSC1 ∪ NSC2, NSC1 ∩ NSC2
(2) PacBio transcripts by type, with/without multigene
(3) Merge iReckon onto PacBio transcripts and count per category, for each of the sets in (1)

*/


/* Get IReckon transcript BLAST results -- best hits:

   Complete hit
   Mostly-complete hit (5-10bp difference in hit length)
   90% hit
   
*/

data blast_hits_100 blast_hits_95 blast_hits_90 blast_hits_75 blast_hits_50 oops;
  set event.ireckon2pacbio_blast_best;
  if flag_len_match_eq_ref=1 or flag_len_match_eq_ref_mm5=1 or flag_len_query_eq_ref=1
      or flag_len_query_eq_ref_mm5=1
     then output blast_hits_100;

  if flag_len_match_5bp_diff_ref=1 or flag_len_match_5bp_diff_ref_mm5=1 
     or flag_len_match_10bp_diff_ref=1 or flag_len_match_10bp_diff_ref_mm5=1
     or flag_len_match_eq_ref=1 or flag_len_match_eq_ref_mm5=1
     or flag_len_query_eq_ref=1 or flag_len_query_eq_ref_mm5=1
     then output blast_hits_95;

  if flag_len_match_5bp_diff_ref=1 or flag_len_match_5bp_diff_ref_mm5=1
     or flag_len_match_10bp_diff_ref=1 or flag_len_match_10bp_diff_ref_mm5=1
     or flag_len_match_eq_ref=1 or flag_len_match_eq_ref_mm5=1
     or flag_len_match_90prc_ref=1 or flag_len_match_90prc_ref_mm5=1
     or flag_len_query_eq_ref=1 or flag_len_query_eq_ref_mm5=1
     then output blast_hits_90;

  if flag_len_match_5bp_diff_ref=1 or flag_len_match_5bp_diff_ref_mm5=1
     or flag_len_match_10bp_diff_ref=1 or flag_len_match_10bp_diff_ref_mm5=1
     or flag_len_match_eq_ref=1 or flag_len_match_eq_ref_mm5=1
     or flag_len_match_90prc_ref=1 or flag_len_match_90prc_ref_mm5=1
     or flag_len_query_eq_ref=1 or flag_len_query_eq_ref_mm5=1
     or flag_len_match_75prc_ref=1 or flag_len_match_75prc_ref_mm5=1
     then output blast_hits_75;

  if flag_len_match_5bp_diff_ref=1 or flag_len_match_5bp_diff_ref_mm5=1
     or flag_len_match_10bp_diff_ref=1 or flag_len_match_10bp_diff_ref_mm5=1
     or flag_len_match_eq_ref=1 or flag_len_match_eq_ref_mm5=1
     or flag_len_match_90prc_ref=1 or flag_len_match_90prc_ref_mm5=1
     or flag_len_query_eq_ref=1 or flag_len_query_eq_ref_mm5=1
     or flag_len_match_50prc_ref=1 or flag_len_match_50prc_ref_mm5=1
     then output blast_hits_50;

  else output oops;

  keep transcript_id query_id sample;
run;

proc sort data=blast_hits_100 nodup;
  by sample query_id transcript_id;
proc sort data=blast_hits_95 nodup;
  by sample query_id transcript_id;
proc sort data=blast_hits_90 nodup;
  by sample query_id transcript_id;
proc sort data=blast_hits_75 nodup;
  by sample query_id transcript_id;
proc sort data=blast_hits_50 nodup;
  by sample query_id transcript_id;
run;

/* Categorize BLAST hits:  NSC1, NSC2, NSC1 ∪ NSC2, NSC1 ∩ NSC2 */

%macro findOverlap(hitsIn,suffixOut);

data ir2pb_nsc1 ir2pb_nsc2;
   set &hitsIn.;
   if sample="NSC1" then output ir2pb_nsc1;
   else if sample="NSC2" then output ir2pb_nsc2;
   keep transcript_id;
run;

proc sort data=ir2pb_nsc1 nodup;
  by transcript_id;
proc sort data=ir2pb_nsc2 nodup;
  by transcript_id;
run;

data pb_hits_overlap_&suffixOut.;
  merge ir2pb_nsc1 (in=in1) ir2pb_nsc2 (in=in2);
  by transcript_id;
  if in1 then flag_in_NPC1_&suffixOut.=1; else flag_in_NPC1_&suffixOut.=0;
  if in2 then flag_in_NPC2_&suffixOut.=1; else flag_in_NPC2_&suffixOut.=0;
  if in1 and in2 then flag_in_NPC1_and_NPC2_&suffixOut.=1; else flag_in_NPC1_and_NPC2_&suffixOut.=0;
  if in1 or in2 then flag_in_NPC1_or_NPC2_&suffixOut.=1; else flag_in_NPC1_or_NPC2_&suffixOut.=0;
run;

%mend;

%findOverlap(blast_hits_100,hit100);
%findOverlap(blast_hits_95,hit95);
%findOverlap(blast_hits_90,hit90);
%findOverlap(blast_hits_75,hit75);
%findOverlap(blast_hits_50,hit50);

proc sort data=pb_hits_overlap_hit100;
  by transcript_id;
proc sort data=pb_hits_overlap_hit95;
  by transcript_id;
proc sort data=pb_hits_overlap_hit90;
  by transcript_id;
proc sort data=pb_hits_overlap_hit75;
  by transcript_id;
proc sort data=pb_hits_overlap_hit50;
  by transcript_id;
run;

data pb_hits_overlap;
  merge pb_hits_overlap_hit100 (in=in1) pb_hits_overlap_hit95 (in=in2)
        pb_hits_overlap_hit90 (in=in3) pb_hits_overlap_hit75 (in=in4)
        pb_hits_overlap_hit50 (in=in5);
  by transcript_id;
run;

data pb_hits_overlap2;
  set pb_hits_overlap;
  array change _numeric_;
  do over change;
     if change=. then change=0;
     end;
run;


/* PacBio transcripts by type, with/without multigene */

data pb_xs;
   set event.pacbio_transcripts;
run;

data pb_xs_nomult;
   set event.pacbio_transcripts_nomulti;
   keep pacbio_id;
run;

proc sort data=pb_xs;
  by pacbio_id;
proc sort data=pb_xs_nomult;
  by pacbio_id;
run;

data pb_xscripts;
  merge pb_xs (in=in1) pb_xs_nomult (in=in2);
  by pacbio_id;
  if in2 then flag_multigene=0; else flag_multigene=1;
run;

/* Merge in BLAST hits */

data pb_hits_overlap_id;
  length pacbio_id $10.;
  set pb_hits_overlap2;
  pacbio_id=compress(scan(transcript_id,1,"|"));
  drop transcript_id;
run;

proc sort data=pb_hits_overlap_id;
   by pacbio_id;
proc sort data=pb_xscripts;
  by pacbio_id;
run;

data pb_xscripts_w_ireckon;
  merge pb_xscripts (in=in1) pb_hits_overlap_id (in=in2);
  by pacbio_id;
run;

data pb_xscripts_w_ireckon2;
  set pb_xscripts_w_ireckon;
  array change _numeric_;
  do over change;
    if change=. then change=0;
    end;
run;

proc freq data=pb_xscripts_w_ireckon2 noprint;
   tables flag_multigene*pb_status*flag_in_NPC1_hit100*flag_in_NPC2_hit100 / out=pb_hits_100perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit95*flag_in_NPC2_hit95 / out=pb_hits_95perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit90*flag_in_NPC2_hit90 / out=pb_hits_90perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit75*flag_in_NPC2_hit75 / out=pb_hits_75perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit50*flag_in_NPC2_hit50 / out=pb_hits_50perc;
run;

/* Export these for counts tables */

proc export data=pb_hits_100perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_100perc_match_summary.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_95perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_95perc_match_summary.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_90perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_90perc_match_summary.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_75perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_75perc_match_summary.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_50perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_50perc_match_summary.csv"
     dbms=csv replace;
run;

/* Now count for only the 5880 we analyzed */


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
   keep pacbio_id splice_match_id match_type;
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
  keep pacbio_id transcript_id;
run;

proc sort data=xs_w_pb nodup;
  by pacbio_id;
proc sort data=pb_xscripts_w_ireckon2;
  by pacbio_id;
run;

data pb_xscripts_w_ireckon3;
  merge pb_xscripts_w_ireckon2 (in=in1) xs_w_pb (in=in2);
  by pacbio_id;
  if in1 and in2;
run;


proc freq data=pb_xscripts_w_ireckon3 noprint;
   tables flag_multigene*pb_status*flag_in_NPC1_hit100*flag_in_NPC2_hit100 / out=pb_hits_100perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit95*flag_in_NPC2_hit95 / out=pb_hits_95perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit90*flag_in_NPC2_hit90 / out=pb_hits_90perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit75*flag_in_NPC2_hit75 / out=pb_hits_75perc;
   tables flag_multigene*pb_status*flag_in_NPC1_hit50*flag_in_NPC2_hit50 / out=pb_hits_50perc;
run;



proc export data=pb_hits_100perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_100perc_match_summary_refseq_match_only.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_95perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_95perc_match_summary_refseq_match_only.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_90perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_90perc_match_summary_refseq_match_only.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_75perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_75perc_match_summary_refseq_match_only.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_50perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_50perc_match_summary_refseq_match_only.csv"
     dbms=csv replace;
run;

proc sort data=pb_xscripts_w_ireckon3 ;
   by transcript_id;
proc means data=pb_xscripts_w_ireckon3 noprint;
   by transcript_id;
   var flag_multigene flag_in_NPC1_hit100 flag_in_NPC2_hit100
       flag_in_NPC1_hit95 flag_in_NPC2_hit95
       flag_in_NPC1_hit90 flag_in_NPC2_hit90
       flag_in_NPC1_hit75 flag_in_NPC2_hit75
       flag_in_NPC1_hit50 flag_in_NPC2_hit50 ;
   output out=pb_xscripts_w_ireckon4 max=;
run;



proc freq data=pb_xscripts_w_ireckon4 noprint;
   tables flag_multigene*flag_in_NPC1_hit100*flag_in_NPC2_hit100 / out=pb_hits_100perc;
   tables flag_multigene*flag_in_NPC1_hit95*flag_in_NPC2_hit95 / out=pb_hits_95perc;
   tables flag_multigene*flag_in_NPC1_hit90*flag_in_NPC2_hit90 / out=pb_hits_90perc;
   tables flag_multigene*flag_in_NPC1_hit75*flag_in_NPC2_hit75 / out=pb_hits_75perc;
   tables flag_multigene*flag_in_NPC1_hit50*flag_in_NPC2_hit50 / out=pb_hits_50perc;
run;



proc export data=pb_hits_100perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_100perc_match_summary_refseq_match_only_collapsed.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_95perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_95perc_match_summary_refseq_match_only_collapsed.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_90perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_90perc_match_summary_refseq_match_only_collapsed.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_75perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_75perc_match_summary_refseq_match_only_collapsed.csv"
     dbms=csv replace;
run;

proc export data=pb_hits_50perc
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_transcripts_with_ireckon_hits_50perc_match_summary_refseq_match_only_collapsed.csv"
     dbms=csv replace;
run;

