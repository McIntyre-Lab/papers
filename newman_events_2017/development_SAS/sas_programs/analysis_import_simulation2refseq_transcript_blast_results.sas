/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;


/* Import BLAST results for simulated transcripts against RefSeq */

%macro importBLAST(sim,test);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sim._&test._transcripts_blast.tsv"
   out=&sim._&test._blast dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data  &sim._&test._blast2;
  length simulation $4.;
  length test $5.;
  set  &sim._&test._blast;
  simulation="&sim.";
  test="&test.";
  rename VAR1=query_id VAR2=transcript_id VAR3=perc_identity VAR4=hit_length VAR5=n_mismatch
         VAR6=n_gapopen VAR7=query_start VAR8=query_stop VAR9=ref_start VAR10=ref_stop
         VAR11=evalue VAR12=bitscore;
run;

%mend;

%importBLAST(sim1,test1);
%importBLAST(sim1,test2);
%importBLAST(sim2,test1);
%importBLAST(sim2,test2);
%importBLAST(sim3,test1);
%importBLAST(sim3,test2);

/* Import transcript sequence lengths */

%macro importLen(sim,test);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sim._&test._transcript_lengths.bed"
   out=&sim._&test._len dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data &sim._&test._len2;
  length simulation $4.;
  length test $5.;
  set &sim._&test._len;
  simulation="&sim.";
  test="&test.";
  drop VAR2;
  rename VAR1=query_id VAR3=query_seq_length;
run;

%mend;

%importLen(sim1,test1);
%importLen(sim1,test2);
%importLen(sim2,test1);
%importLen(sim2,test2);
%importLen(sim3,test1);
%importLen(sim3,test2);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/refseq_transcripts_length.bed"
   out=refseq_len dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data refseq_len2;
  set refseq_len;
  drop VAR2;
  rename VAR1=transcript_id VAR3=transcript_seq_length;
run;

/* Stack, merge and make permenant */

data all_blast;
  set sim1_test1_blast2 sim1_test2_blast2
      sim2_test1_blast2 sim2_test2_blast2
      sim3_test1_blast2 sim3_test2_blast2;
run;

data all_lengths;
  set sim1_test1_len2 sim1_test2_len2
      sim2_test1_len2 sim2_test2_len2
      sim3_test1_len2 sim3_test2_len2;
run;

proc sort data=all_blast;
  by simulation test query_id;
proc sort data=all_lengths;
  by simulation test query_id;
run;

data all_blast_w_len;
  merge all_blast (in=in1) all_lengths (in=in2);
  by simulation test query_id;
  if in1 and in2;
run;

proc sort data=all_blast_w_len;
  by transcript_id;
proc sort data=refseq_len2;
  by transcript_id;
run;

data all_blast_w_len2;
  merge all_blast_w_len (in=in1) refseq_len2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Make permenant */

data event.sim_xscript2refseq_blast;
  set all_blast_w_len2;
run;



