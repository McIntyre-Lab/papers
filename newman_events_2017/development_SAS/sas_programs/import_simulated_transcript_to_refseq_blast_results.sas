/* Libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import BLAST results for simulated data to RefSeq */

%macro importBLAST(sim,test);

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/&sim._&test._transcripts_blast.tsv"
     out=&sim._&test._blast dbms=tab replace;
     guessingrows=max; getnames=no;
run;

data &sim._&test._blast2;
  length simulation $4.;
  length test $5.;
  set &sim._&test._blast;
  simulation="&sim.";
  test="&test.";
  rename VAR1=query_id
         VAR2=transcript_id
         VAR3=perc_identity
         VAR4=length
         VAR5=num_mismatch
         VAR6=num_gapopen
         VAR7=query_start
         VAR8=query_stop
         VAR9=ref_start
         VAR10=ref_stop
         VAR11=evalue
         VAR12=bitscore;
run;

%mend;

%importBLAST(sim1,test1);
%importBLAST(sim2,test1);
%importBLAST(sim3,test1);
%importBLAST(sim1,test2);
%importBLAST(sim2,test2);
%importBLAST(sim3,test2);

/* Stack results */

data sim_xs2ref;
   set sim1_test1_blast2 sim2_test1_blast2 sim3_test1_blast2
       sim1_test2_blast2 sim2_test2_blast2 sim3_test2_blast2;
run;


/* Make permenant */

data event.simulated_xscript2refseq_blast;
  set sim_xs2ref;
run;

