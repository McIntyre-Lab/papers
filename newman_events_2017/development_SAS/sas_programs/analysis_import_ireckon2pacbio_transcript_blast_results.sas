/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;


/* Import BLAST results for iReckon transcripts against PacBio */

%macro importBLAST(sample);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sample._ireckon_to_pacbio_blast_hits.tsv"
   out=&sample._blast dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data  &sample._blast2;
  length sample $4.;
  set  &sample._blast;
  sample="&sample.";
  rename VAR1=query_id VAR2=transcript_id VAR3=perc_identity VAR4=hit_length VAR5=n_mismatch
         VAR6=n_gapopen VAR7=query_start VAR8=query_stop VAR9=ref_start VAR10=ref_stop
         VAR11=evalue VAR12=bitscore;
run;

%mend;

%importBLAST(NSC1);
%importBLAST(NSC2);

/* Import transcript sequence lengths */

%macro importLen(sample);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sample._iReckon_transcripts_length.bed"
   out=&sample._len dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data &sample._len2;
  length sample $4.;
  set &sample._len;
  sample="&sample.";
  drop VAR2;
  rename VAR1=query_id VAR3=query_seq_length;
run;

%mend;

%importLen(NSC1);
%importLen(NSC2);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/pacbio_transcripts_length.bed"
   out=pacbio_len dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data pacbio_len2;
  set pacbio_len;
  drop VAR2;
  rename VAR1=transcript_id VAR3=transcript_seq_length;
run;

/* Stack, merge and make permenant */

data all_blast;
  length query_id $53.;
  set NSC1_blast2 NSC2_blast2;
run;

data all_lengths;
  length query_id $53.;
  set NSC1_len2 NSC2_len2;
run;

proc sort data=all_blast;
  by sample query_id;
proc sort data=all_lengths;
  by sample query_id;
run;

data all_blast_w_len;
  length transcript_id $40.;
  format transcript_id $40.;
  merge all_blast (in=in1) all_lengths (in=in2);
  by sample query_id;
  if in1 and in2;
run;

proc sort data=all_blast_w_len;
  by transcript_id;
proc sort data=pacbio_len2;
  by transcript_id;
run;

data all_blast_w_len2;
  merge all_blast_w_len (in=in1) pacbio_len2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Make permenant */

data event.ireckon2pacbio_blast;
  set all_blast_w_len2;
run;



