/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Import junctions from STAR alignment output for each of the simulated datasets to count if the junctions it detects
   are in the junction catalog */

%macro importSJ(sample);

proc import datafile=""
      out=&sample._junc dbms=tab replace;
      guessingrows=176000; getnames=no;
run;

data &sample._junc2;
  length sample_id $15.;
  set &sample._junc;
  sample_id="&sample.";
  rename VAR1=chr
         VAR2=intron_start
         VAR3=intron_stop
         VAR4=strand
         VAR5=intron_motif_type
         VAR6=flag_junction_annotated
         VAR7=num_unique_mapped_reads
         VAR8=num_multimapped_reads
         VAR9=max_overhang
         ;
run;

%mend;

%importSJ(sim1_test1);
%importSJ(sim1_test2);
%importSJ(sim2_test1);
%importSJ(sim2_test2);
%importSJ(sim3_test1);
%importSJ(sim3_test2);

/* Stack junctions */

data star_junc;
  set sim1_test1_junc2 sim1_test2_junc2 sim2_test1_junc2 sim2_test2_junc2 sim3_test1_junc2 sim3_test2_junc2;
run;

/* Put unique-mapped counts side-by-side */

proc sort data=star_junc;
  by chr intron_Start intron_stop strand intron_motif_type flag_junction_annotated sample_id;
proc transpose data=star_junc out=star_junc_sbys(drop=_NAME_);
  by chr intron_Start intron_stop strand intron_motif_type flag_junction_annotated ;
  var num_unique_mapped_reads;
  id sample_id;
run;

data star_junc_sbys2;
  set star_junc_sbys;
  array change _numeric_;
     do over change;
     if change=. then change=0;
     end;
run;

/* Flag if junction has at least 5 reads in a given sample */

data flag_star_junc;
  retain chr intron_start intron_stop strand intron_motif_type flag_junction_annotated
         sim1_test1 sim2_test1 sim3_test1 sim1_test2 sim2_test2 sim3_test2;
  set star_junc_sbys2;
  if sim1_test1 ge 5 then flag_depth_sim1_test1_ge5=1; else flag_depth_sim1_test1_ge5=0;
  if sim2_test1 ge 5 then flag_depth_sim2_test1_ge5=1; else flag_depth_sim2_test1_ge5=0;
  if sim3_test1 ge 5 then flag_depth_sim3_test1_ge5=1; else flag_depth_sim3_test1_ge5=0;
  if sim1_test2 ge 5 then flag_depth_sim1_test2_ge5=1; else flag_depth_sim1_test2_ge5=0;
  if sim2_test2 ge 5 then flag_depth_sim2_test2_ge5=1; else flag_depth_sim2_test2_ge5=0;
  if sim3_test2 ge 5 then flag_depth_sim3_test2_ge5=1; else flag_depth_sim3_test2_ge5=0;
  if flag_depth_sim1_test1_ge5=1
     and flag_depth_sim2_test1_ge5=1
     and flag_depth_sim3_test1_ge5=1 then flag_test1_all_ge5=1; else flag_test1_all_ge5=0;

  if flag_depth_sim1_test2_ge5=1
     and flag_depth_sim2_test2_ge5=1
     and flag_depth_sim3_test2_ge5=1 then flag_test2_all_ge5=1; else flag_test2_all_ge5=0;
run;

/* Count : how many junctions detected at depth > 5? How mant detected in all? */

proc freq data=flag_star_junc;
  tables flag_depth_sim1_test1_ge5 flag_depth_sim2_test1_ge5 flag_depth_sim3_test1_ge5 flag_test1_all_ge5
         flag_depth_sim1_test2_ge5 flag_depth_sim2_test2_ge5 flag_depth_sim3_test2_ge5 flag_test2_all_ge5;
run;

/*
 flag_depth_sim1_                             Cumulative    Cumulative
        test1_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      378794       79.00        378794        79.00
                1      100665       21.00        479459       100.00


 flag_depth_sim2_                             Cumulative    Cumulative
        test1_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      376538       78.53        376538        78.53
                1      102921       21.47        479459       100.00


 flag_depth_sim3_                             Cumulative    Cumulative
        test1_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      378446       78.93        378446        78.93
                1      101013       21.07        479459       100.00

    flag_test1_                             Cumulative    Cumulative
        all_ge5    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0      445185       92.85        445185        92.85
              1       34274        7.15        479459       100.00


lag_depth_sim1_                             Cumulative    Cumulative
      test2_ge5    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0      377230       78.68        377230        78.68
              1      102229       21.32        479459       100.00


lag_depth_sim2_                             Cumulative    Cumulative
      test2_ge5    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0      377817       78.80        377817        78.80
              1      101642       21.20        479459       100.00

 flag_depth_sim3_                             Cumulative    Cumulative
        test2_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      374099       78.03        374099        78.03
                1      105360       21.97        479459       100.00


      flag_test2_                             Cumulative    Cumulative
          all_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      447788       93.39        447788        93.39
                1       31671        6.61        479459       100.00
*/
   

data event.simul_star_junctions;
  set flag_star_junc;
run;

