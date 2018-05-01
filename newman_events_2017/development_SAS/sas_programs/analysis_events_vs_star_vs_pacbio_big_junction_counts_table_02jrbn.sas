/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;


/* Create counts table by junction type. Want the following columns:
  Junction types: annotated, NIC, NNC novel donor/acceptor, NNC novel donor, NNC novel acceptor, NNC other
  PacBio counts
  STAR "detected" counts
  Junction catalog total
  Junction catalog "detected" counts
  STAR support (APN0,5,10) -> W/WO
  Catalog support (APN0,5,10) -> W/WO

  I will be exporting these to make 2 tables: one of PacBio junctions, one of non-PB junctions

  I also want to check the concordance in APN between STAR and aligning to the catalog, but this is another step
*/

/* (1) redo apn flags for junction catalog */

data flag_apn;
  set eventloc.unique_junction2event_mm10;
  if NSC1_apn > 0 or NSC2_apn > 0 then flag_events_detected=1; else flag_events_detected=0;
  if NSC1_apn > 0 then flag_events_NSC1_apn_gt0=1; else flag_events_NSC1_apn_gt0=0;
  if NSC2_apn > 0 then flag_events_NSC2_apn_gt0=1; else flag_events_NSC2_apn_gt0=0;

  if NSC1_apn ge 2 then flag_events_NSC1_apn_ge2=1; else flag_events_NSC1_apn_ge2=0;
  if NSC2_apn ge 2 then flag_events_NSC2_apn_ge2=1; else flag_events_NSC2_apn_ge2=0;

  if NSC1_apn ge 5 then flag_events_NSC1_apn_ge5=1; else flag_events_NSC1_apn_ge5=0;
  if NSC2_apn ge 5 then flag_events_NSC2_apn_ge5=1; else flag_events_NSC2_apn_ge5=0;

  if NSC1_apn ge 10 then flag_events_NSC1_apn_ge10=1; else flag_events_NSC1_apn_ge10=0;
  if NSC2_apn ge 10 then flag_events_NSC2_apn_ge10=1; else flag_events_NSC2_apn_ge10=0;

  if flag_events_NSC1_apn_gt0=1 and flag_events_NSC2_apn_gt0=1 then flag_events_NSC_apn_gt0=1;
  else if flag_events_NSC1_apn_gt0=1 or flag_events_NSC2_apn_gt0=1 then flag_events_NSC_apn_gt0=.;
  else flag_events_NSC_apn_gt0=0;

  if flag_events_NSC1_apn_ge2=1 and flag_events_NSC2_apn_ge2=1 then flag_events_NSC_apn_ge2=1;
  else if flag_events_NSC1_apn_ge2=1 or flag_events_NSC2_apn_ge2=1 then flag_events_NSC_apn_ge2=.;
  else flag_events_NSC_apn_ge2=0;

  if flag_events_NSC1_apn_ge5=1 and flag_events_NSC2_apn_ge5=1 then flag_events_NSC_apn_ge5=1;
  else if flag_events_NSC1_apn_ge5=1 or flag_events_NSC2_apn_ge5=1 then flag_events_NSC_apn_ge5=.;
  else flag_events_NSC_apn_ge5=0;

  if flag_events_NSC1_apn_ge10=1 and flag_events_NSC2_apn_ge10=1 then flag_events_NSC_apn_ge10=1;
  else if flag_events_NSC1_apn_ge10=1 or flag_events_NSC2_apn_ge10=1 then flag_events_NSC_apn_ge10=.;
  else flag_events_NSC_apn_ge10=0;

  rename NSC1_apn=NSC1_apn_events NSC2_apn=NSC2_apn_events;
run;

/* (2) Get STAR APN counts and flag */

data apn_star;
   set event.npc_star_junctions_est_apn;
   keep chr intron_start intron_stop strand NPC1_apn NPC2_apn;
run;

data flag_apn_star;
  set apn_star;
  if NPC1_apn > 0 or NPC2_apn > 0 then flag_star_detected=1; else flag_star_detected=0;
  if NPC1_apn > 0 then flag_star_NSC1_apn_gt0=1; else flag_star_NSC1_apn_gt0=0;
  if NPC2_apn > 0 then flag_star_NSC2_apn_gt0=1; else flag_star_NSC2_apn_gt0=0;

  if NPC1_apn ge 2 then flag_star_NSC1_apn_ge2=1; else flag_star_NSC1_apn_ge2=0;
  if NPC2_apn ge 2 then flag_star_NSC2_apn_ge2=1; else flag_star_NSC2_apn_ge2=0;

  if NPC1_apn ge 5 then flag_star_NSC1_apn_ge5=1; else flag_star_NSC1_apn_ge5=0;
  if NPC2_apn ge 5 then flag_star_NSC2_apn_ge5=1; else flag_star_NSC2_apn_ge5=0;

  if NPC1_apn ge 10 then flag_star_NSC1_apn_ge10=1; else flag_star_NSC1_apn_ge10=0;
  if NPC2_apn ge 10 then flag_star_NSC2_apn_ge10=1; else flag_star_NSC2_apn_ge10=0;

  if flag_star_NSC1_apn_gt0=1 and flag_star_NSC2_apn_gt0=1 then flag_star_NSC_apn_gt0=1;
  else if flag_star_NSC1_apn_gt0=1 or flag_star_NSC2_apn_gt0=1 then flag_star_NSC_apn_gt0=.;
  else flag_star_NSC_apn_gt0=0;

  if flag_star_NSC1_apn_ge2=1 and flag_star_NSC2_apn_ge2=1 then flag_star_NSC_apn_ge2=1;
  else if flag_star_NSC1_apn_ge2=1 or flag_star_NSC2_apn_ge2=1 then flag_star_NSC_apn_ge2=.;
  else flag_star_NSC_apn_ge2=0;

  if flag_star_NSC1_apn_ge5=1 and flag_star_NSC2_apn_ge5=1 then flag_star_NSC_apn_ge5=1;
  else if flag_star_NSC1_apn_ge5=1 or flag_star_NSC2_apn_ge5=1 then flag_star_NSC_apn_ge5=.;
  else flag_star_NSC_apn_ge5=0;

  if flag_star_NSC1_apn_ge10=1 and flag_star_NSC2_apn_ge10=1 then flag_star_NSC_apn_ge10=1;
  else if flag_star_NSC1_apn_ge10=1 or flag_star_NSC2_apn_ge10=1 then flag_star_NSC_apn_ge10=.;
  else flag_star_NSC_apn_ge10=0;
  donor_stop=intron_start-1;
  acceptor_start=intron_stop;
   if strand=2 then strand_chr="-"; else strand_chr="+";
   drop strand;
  rename NPC1_apn=NSC1_apn_star NPC2_apn=NSC2_apn_star  strand_chr=strand;
run;

/* Merge */

data pb_junc;
   length chr $20.;
   set evspl.splicing_events_annotations;
   where num_transcripts > 0;
   keep chr strand feature1_stop feature2_start;
   rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;

data pb_junc2;
  set event.catalog_pacbio_star_junctions;
  where flag_in_pacbio=1;
  keep chr strand donor_Stop acceptor_start
       flag_donor_in_catalog flag_acceptor_in_catalog;
run;

proc sort data=pb_junc nodup;
   by chr strand donor_stop acceptor_start;
proc sort data=pb_junc2;
   by chr strand donor_stop acceptor_start;
run;

data pb_junc3;
  merge pb_junc (in=in1) pb_junc2 (in=in2);
  by chr strand donor_stop acceptor_start;
  if in1 and in2;
run;

data cat_junc;
   length chr $20.;
  set evspl.splicing_events_annot_refseq;
  if num_transcripts>0 then flag_junction_annotated=1;
  keep event_id chr strand feature1_stop feature2_start flag_intron_retention flag_junction_annotated;
  rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;

proc sort data=cat_junc;
   by event_id;
proc sort data=flag_apn;
   by event_id;
run;

data cat_junc_w_apn;
  merge cat_junc (in=in1) flag_apn (in=in2);
  by event_id;
  if in1 and in2;
run;


proc sort data=cat_junc_w_apn nodup;
   by chr strand donor_stop acceptor_start;
run;
proc sort data=pb_junc3 nodup;
   by chr strand donor_stop acceptor_start;
run;
proc sort data=flag_apn_star nodup;
   by chr strand donor_stop acceptor_start;
run;

data junc_table;
  merge cat_junc_w_apn (in=in1) pb_junc3 (in=in2) flag_apn_star (in=in3);
  by chr strand donor_stop acceptor_start;
  if in1 then flag_in_catalog=1; else flag_in_catalog=0;
  if in2 then flag_in_pacbio=1; else flag_in_pacbio=0; 
  if in3 then flag_in_star=1; else flag_in_star=0;  
run;

data junc_table2;
   length junction_type $50.;
   set junc_table;
   if flag_intron_retention=1 then junction_type="border junction";
   else do;
      if flag_junction_annotated=1 then junction_type="Annotated";
      else if flag_junction_annotated=0 then junction_type="NIC/Unannotated";
      else if flag_junction_annotated=. then do;
         if flag_donor_in_catalog=1 and flag_acceptor_in_catalog=1 then junction_type="NNC, known donor and acceptor";
         else if flag_donor_in_catalog=1 and flag_acceptor_in_catalog=0 then junction_type="NNC, novel acceptor";
         else if flag_donor_in_catalog=0 and flag_acceptor_in_catalog=1 then junction_type="NNC, novel donor";
         else junction_type="NNC, novel donor and acceptor";
      end;
   end;
   keep chr event_id strand Donor_stop acceptor_Start flag_junction_annotated junction_type
   flag_intron_retention flag_donor_in_catalog flag_acceptor_in_catalog flag_in_catalog
   flag_in_pacbio flag_in_star;
run;


proc freq data=junc_table2;
   tables junction_type;
run;

/*

                                                            Cumulative    Cumulative
  junction_type                    Frequency     Percent     Frequency      Percent
  ----------------------------------------------------------------------------------
  Annotated                          290055        9.61        290055         9.61
  NIC/Unannotated                   2460416       81.52       2750471        91.14
  NNC, known donor and acceptor          38        0.00       2750509        91.14
  NNC, novel acceptor                   928        0.03       2751437        91.17
  NNC, novel donor                      962        0.03       2752399        91.20
  NNC, novel donor and acceptor       17631        0.58       2770030        91.78
  border junction                    247983        8.22       3018013       100.00

*/

proc freq data=junc_table2 noprint;
   tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star / out=junc_by_type_tech;
proc print data=junc_by_type_tech;
run;

/*
                                 flag_in_    flag_in_    flag_in_
junction_type                     pacbio      catalog      star        COUNT

Annotated                            0           1           0        173772
Annotated                            0           1           1         47890
Annotated                            1           1           0          7247
Annotated                            1           1           1         61146
NIC/Unannotated                      0           1           0       2456029
NIC/Unannotated                      0           1           1          3735
NIC/Unannotated                      1           1           0           333
NIC/Unannotated                      1           1           1           319
NNC, known donor and acceptor        1           0           0            23
NNC, known donor and acceptor        1           0           1            15
NNC, novel acceptor                  1           0           0           667
NNC, novel acceptor                  1           0           1           261
NNC, novel donor                     1           0           0           703
NNC, novel donor                     1           0           1           259
NNC, novel donor and acceptor        0           0           1         15868
NNC, novel donor and acceptor        1           0           0          1644
NNC, novel donor and acceptor        1           0           1           119
border junction                      0           1           0        247983


*/

proc sort data=junc_table2;
   by event_id;
proc sort data=flag_apn;
   by event_id;
run;

data junc_tables_event_apn;
   merge junc_table2 (in=in1) flag_apn (in=in2);
   by event_id;
   if in1;
run;

proc sort data=junc_tables_event_apn;
   by chr strand donor_stop acceptor_start;
proc sort data=flag_apn_star;
   by chr strand donor_stop acceptor_start;
run;

data junc_tables_star_apn;
   merge junc_tables_event_apn (in=in1) flag_apn_star (in=in2);
   by chr strand donor_stop acceptor_start;
   if flag_star_detected=. then flag_star_detected=0;
   if flag_events_detected=. then flag_events_detected=0;
   if in1 then output;
run;




/* Create individual tables for APN0, APN2, APN5, APN10 */

proc freq data=junc_tables_star_apn noprint;
  tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star*flag_star_detected*flag_events_detected 
         / out=junc_count_dtct;

  tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star*flag_star_detected*flag_events_detected*
         flag_events_NSC_apn_gt0*flag_star_NSC_apn_gt0 / out=junc_count_apn0;

  tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star*flag_star_detected*flag_events_detected*
         flag_events_NSC_apn_ge2*flag_star_NSC_apn_ge2 / out=junc_count_apn2;

  tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star*flag_star_detected*flag_events_detected*
         flag_events_NSC_apn_ge5*flag_star_NSC_apn_ge5 / out=junc_count_apn5;

  tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star*flag_star_detected*flag_events_detected*
         flag_events_NSC_apn_ge10*flag_star_NSC_apn_ge10 / out=junc_count_apn10;

run;

/* Export */

proc export data=junc_count_dtct
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_detection.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn0
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn0.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn2
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn2.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn5
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn5.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn10
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn10.csv"
     dbms=csv replace;
run;


/* Count junctions by type in PB, catalog, STAR */

data junc_tables_star_apn2;
   set junc_tables_star_apn;
   array change _numeric_;
   do over change;
      if change=. then change=0;
   end;
run;

/* Make permenant */

data event.event2star2pacbio_junc_table;
  set junc_tables_star_apn2;
run;


/* Totals -- overall */

proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_in_catalog
          junction_type*flag_in_star;
run;

/* Totals -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_in_catalog
          junction_type*flag_in_star;
run;

/* Totals -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_in_catalog=1;
   tables junction_type*flag_in_star;
run;

/* total -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_in_catalog=1;
   tables junction_type*flag_in_star;
run;




/* Detection -- overall */

proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_events_detected 
          junction_type*flag_star_detected;
run;

/* Detection -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_events_detected
          junction_type*flag_star_detected;
run;

/* Detection -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_events_detected=1;
   tables junction_type*flag_star_detected;
run;

/* Detection -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_events_detected=1;
   tables junction_type*flag_star_detected;
run;



/* APN0 -- overall */
proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_events_NSC_apn_gt0 
          junction_type*flag_star_NSC_apn_gt0;
run;

/* Detection -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_events_NSC_apn_gt0
          junction_type*flag_star_NSC_apn_gt0;
run;

/* Detection -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_events_NSC_apn_gt0=1;
   tables junction_type*flag_star_NSC_apn_gt0;
run;

/* Detection -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_events_NSC_apn_gt0=1;
   tables junction_type*flag_star_NSC_apn_gt0;
run;


/* APN2 -- overall */

proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_events_NSC_apn_ge2 
          junction_type*flag_star_NSC_apn_ge2;
run;

/* Detection -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_events_NSC_apn_ge2
          junction_type*flag_star_NSC_apn_ge2;
run;

/* Detection -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_events_NSC_apn_ge2=1;
   tables junction_type*flag_star_NSC_apn_ge2;
run;

/* Detection -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_events_NSC_apn_ge2=1;
   tables junction_type*flag_star_NSC_apn_ge2;
run;

/* APN5 -- overall */

proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_events_NSC_apn_ge5 
          junction_type*flag_star_NSC_apn_ge5;
run;

/* Detection -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_events_NSC_apn_ge5
          junction_type*flag_star_NSC_apn_ge5;
run;

/* Detection -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_events_NSC_apn_ge5=1;
   tables junction_type*flag_star_NSC_apn_ge5;
run;

/* Detection -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_events_NSC_apn_ge5=1;
   tables junction_type*flag_star_NSC_apn_ge5;
run;

/* APN10 -- overall */

proc freq data=junc_tables_star_apn2;
   tables junction_type*flag_in_pacbio
          junction_type*flag_events_NSC_apn_ge10 
          junction_type*flag_star_NSC_apn_ge10;
run;

/* Detection -- overlap with PacBio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1;
   tables junction_type*flag_events_NSC_apn_ge10
          junction_type*flag_star_NSC_apn_ge10; 
run;

/* Detection -- overlap with catalog */
proc freq data=junc_tables_star_apn2;
   where flag_events_NSC_apn_ge10=1;
   tables junction_type*flag_star_NSC_apn_ge10;
run;

/* Detection -- overlap with catalog and pacbio */
proc freq data=junc_tables_star_apn2;
   where flag_in_pacbio=1 and flag_events_NSC_apn_ge10=1;
   tables junction_type*flag_star_NSC_apn_ge10;
run;



