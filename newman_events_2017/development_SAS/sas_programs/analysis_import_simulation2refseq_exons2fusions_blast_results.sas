/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import simulation exons to Refseq fusions BLAST results */

%macro importBLAST(sim,test);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sim._&test._exons_mm9_to_mm10.tsv"
   out=&sim._&test._blast dbms=tab replace;
   getnames=no; guessingrows=max;
run;

data  &sim._&test._blast2;
  length simulation $4.;
  length test $5.;
  set  &sim._&test._blast;
  simulation="&sim.";
  test="&test.";
  rename VAR1=exon_id VAR2=fusion_id VAR3=perc_identity VAR4=hit_length VAR5=n_mismatch
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
