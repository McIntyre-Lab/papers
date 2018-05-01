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

/* Check if any junction that has a minimum overhang < 16bp */

data flag_overhang;
  set star_junc;
  if max_overhang < 16 then flag_overhang_lt16=1;
  else flag_overhang_lt16=0;
run;

proc freq data=flag_overhang;
   tables flag_overhang_lt16;
run;

data check;
   set flag_overhang;
   where flag_overhang_lt16=1;
run;
/* 107 junctions with invalid overhang. Drop these */

data star_junc2;
  set flag_overhang;
  where flaG_overhang_lt16=0;
run;


/* Put unique-mapped counts side-by-side */

proc sort data=star_junc2;
  by chr intron_Start intron_stop strand intron_motif_type flag_junction_annotated sample_id;
proc transpose data=star_junc2 out=star_junc_sbys(drop=_NAME_);
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
                0       75064       58.80         75064        58.80
                1       52598       41.20        127662       100.00


      flag_depth_                             Cumulative    Cumulative
         NSC2_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       84116       65.89         84116        65.89
                1       43546       34.11        127662       100.00


                                              Cumulative    Cumulative
 flag_NSC_all_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       87002       68.15         87002        68.15
                1       40660       31.85        127662       100.00


      flag_depth_                             Cumulative    Cumulative
         NSC1_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       82911       57.08         82911        57.08
                1       62345       42.92        145256       100.00


      flag_depth_                             Cumulative    Cumulative
         NSC2_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       92888       63.95         92888        63.95
                1       52368       36.05        145256       100.00


                                              Cumulative    Cumulative
 flag_NSC_all_ge5    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       95982       66.08         95982        66.08
                1       49274       33.92        145256       100.00



*/
   

data event.NPC_star_junctions;
  set flag_star_junc;
run;

