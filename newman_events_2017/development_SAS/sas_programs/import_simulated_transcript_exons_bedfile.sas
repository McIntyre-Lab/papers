/* Libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* For each simulated dataset, import the generated exon BED file.

Match exons to Refseq exons (and transcripts). Then for each combination of RefSeq transcript and simulation transcript
count the number of matching exons
*/

%macro importBED(sim,test);

proc import datafile="!MCLAB/event_analysis/analysis_output/&sim._&test._exons.bed"
     out=&sim._&test._exons dbms=tab replace;
     getnames=no; guessingrows=max;
run;

data &sim._&test._exons2;
  length simulation $4.;
  length test $5.;
  length sim_transcript_id $20.;
  set &sim._&test._exons;
  simulation="&sim.";
  test="&test.";
  sim_transcript_id=scan(VAR4,1,":");
  drop VAR5;
  rename VAR1=chr VAR2=exon_start VAR3=exon_stop VAR4=sim_exon_id VAR6=strand;
run;

%mend;

%importBED(sim1,test1);
%importBED(sim1,test2);
%importBED(sim2,test1);
%importBED(sim2,test2);
%importBED(sim3,test1);
%importBED(sim3,test2);

data sim_exons;
   set sim1_test1_exons2 sim2_test1_exons2 sim3_test1_exons2
       sim2_test1_exons2 sim2_test2_exons2 sim3_test2_exons2;
run;

