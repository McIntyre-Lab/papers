/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import simulated exon coordinates */

%macro importBED(sim,test);

     data WORK.&sim._&test._mm9    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "!MCLAB/event_analysis/analysis_output/&sim._&test._exons_mm9.bed" delimiter='09'x
 MISSOVER DSD lrecl=32767 ;
        informat mm9_chr $12. ;
        informat mm9_exon_start best32. ;
        informat mm9_exon_stop best32. ;
        informat sim_exon_id $17. ;
        informat score best32. ;
        informat mm9_strand $1. ;
        format mm9_chr $12. ;
        format mm9_exon_start best12. ;
        format mm9_exon_stop best12. ;
        format sim_exon_id $17. ;
        format score best12. ;
        format mm9_strand $1. ;
     input
                 mm9_chr $
                 mm9_exon_start
                 mm9_exon_stop
                 sim_exon_id $
                 score
                 mm9_strand $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


     data WORK.&sim._&test._mm10    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "!MCLAB/event_analysis/analysis_output/&sim._&test._exons_mm9_to_mm10.bed" delimiter='09'x
 MISSOVER DSD lrecl=32767 ;
        informat chr $12. ;
        informat exon_start best32. ;
        informat exon_stop best32. ;
        informat sim_exon_id $17. ;
        informat score best32. ;
        informat strand $1. ;
        format chr $12. ;
        format exon_start best12. ;
        format exon_stop best12. ;
        format sim_exon_id $17. ;
        format score best12. ;
        format strand $1. ;
     input
                 chr $
                 exon_start
                 exon_stop
                 sim_exon_id $
                 score
                 strand $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data &sim._&test._mm9_2;
  length simulation $4.;
  length test $5.;
  set &sim._&test._mm9;
  simulation="&sim.";
  test="&test.";
  drop score;
run;

data &sim._&test._mm10_2;
  length simulation $4.;
  length test $5.;
  set &sim._&test._mm10;
  simulation="&sim.";
  test="&test.";
  drop score;
run;

proc sort data=&sim._&test._mm9_2;
   by simulation test sim_exon_id;
proc sort data=&sim._&test._mm10_2;
   by simulation test sim_exon_id;
run;

data &sim._&test._exons ;
   merge &sim._&test._mm9_2 (in=in1) &sim._&test._mm10_2 (in=in2);
   by simulation test sim_exon_id;
   if in1 and in2 then flag_mm10_coord=1;
   else if in1 then flag_mm10_coord=0;
   else if in2 then flag_mm10_coord=1;
run;

%mend;

%importBED(sim1,test1);
%importBED(sim1,test2);
%importBED(sim2,test1);
%importBED(sim2,test2);
%importBED(sim3,test1);
%importBED(sim3,test2);

/* Stack and make permenant */

data event.simulated_exons_mm9_to_mm10;
  set sim1_test1_exons sim1_test2_exons
      sim2_test1_exons sim2_test2_exons
      sim3_test1_exons sim3_test2_exons;
run;


