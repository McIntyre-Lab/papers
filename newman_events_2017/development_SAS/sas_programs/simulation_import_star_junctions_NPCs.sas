/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Import junctions from STAR alignment output for each of the simulated datasets to count if the junctions it detects
   are in the junction catalog */

%macro importSJ(sample);

proc import datafile="!MCLAB/event_analysis/alignment_output/star_junctions/&sample.SJ.out.tab"
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

%importSJ(NSC1);
%importSJ(NSC2);

/* Stack junctions */

data star_junc;
  set NSC1_junc2 NSC2_junc2;
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
         NSC1 NSC2;
  set star_junc_sbys2;
  if NSC1 ge 5 then flag_depth_NSC1_ge5=1; else flag_depth_NSC1_ge5=0;
  if NSC2 ge 5 then flag_depth_NSC2_ge5=1; else flag_depth_NSC2_ge5=0;

  if flag_depth_NSC1_ge5=1 and flag_depth_NSC2_ge5=1
  then flag_NSC_all_ge5=1; else flag_NSC_all_ge5=0;
run;

/* Count : how many junctions detected at depth > 5? How mant detected in all? */

proc freq data=flag_star_junc;
  tables flag_depth_NSC1_ge5 flag_depth_NSC2_ge5 flag_NSC_all_ge5;
run;

/*
      flag_depth_                             Cumulative    Cumulative
         NSC1_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      114197       58.54        114197        58.54
                1       80870       41.46        195067       100.00


      flag_depth_                             Cumulative    Cumulative
         NSC2_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      124749       63.95        124749        63.95
                1       70318       36.05        195067       100.00


                                              Cumulative    Cumulative
 flag_NSC_all_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      128333       65.79        128333        65.79
                1       66734       34.21        195067       100.00
*/
   

data event.NPC_star_junctions;
  set flag_star_junc;
run;

