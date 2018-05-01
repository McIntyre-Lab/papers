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

  Limit this to ONLY those genes that we used for event analysis
*/

/* Get list of donor and acceptor sites (coordinates) to delete from big table */

data event2gene;
  set evspl.splicing_events_annot_refseq;
  keep gene_id chr feature1_stop feature2_start strand;
  rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;

data genes2keep;
  set event.feature2xs2gene_exp_only_nomulti;
  keep gene_id;
run;

proc sort data=event2gene;
  by gene_id;
proc sort data=genes2keep nodup;
   by gene_id;
run;
 
data events2drop;
   merge genes2keep (in=in1) event2gene (in=in2);
   by gene_id;
   if in1 and in2 then delete;
   else if in2 then output;
run;


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

data events2drop2;
  set events2drop;
  drop gene_id;
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
proc sort data=events2drop2 nodup;
   by chr strand donor_stop acceptor_start;
run;


data junc_table;
  merge cat_junc_w_apn (in=in1) pb_junc3 (in=in2) flag_apn_star (in=in3) events2drop2 (in=in4);
  by chr strand donor_stop acceptor_start;
  if in1 then flag_in_catalog=1; else flag_in_catalog=0;
  if in2 then flag_in_pacbio=1; else flag_in_pacbio=0; 
  if in3 then flag_in_star=1; else flag_in_star=0;  
  if in4 then delete;
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
  Annotated                          181860        9.34        181860         9.34
  NIC/Unannotated                   1557386       79.96       1739246        89.30
  NNC, known donor and acceptor          38        0.00       1739284        89.30
  NNC, novel acceptor                   928        0.05       1740212        89.35
  NNC, novel donor                      962        0.05       1741174        89.40
  NNC, novel donor and acceptor       51437        2.64       1792611        92.04
  border junction                    154987        7.96       1947598       100.00

*/

proc freq data=junc_table2 noprint;
   tables junction_type*flag_in_pacbio*flag_in_catalog*flag_in_star / out=junc_by_type_tech;
proc print data=junc_by_type_tech;
run;

/*
                                   flag_in_    flag_in_    flag_in_
  junction_type                     pacbio      catalog      star        COUNT

  Annotated                            0           1           0         85089
  Annotated                            0           1           1         49215
  Annotated                            1           1           0          1642
  Annotated                            1           1           1         45914
  NIC/Unannotated                      0           1           0       1550703
  NIC/Unannotated                      0           1           1          6224
  NIC/Unannotated                      1           1           0           180
  NIC/Unannotated                      1           1           1           279
  NNC, known donor and acceptor        1           0           0            15
  NNC, known donor and acceptor        1           0           1            23
  NNC, novel acceptor                  1           0           0           573
  NNC, novel acceptor                  1           0           1           355
  NNC, novel donor                     1           0           0           582
  NNC, novel donor                     1           0           1           380
  NNC, novel donor and acceptor        0           0           1         49674
  NNC, novel donor and acceptor        1           0           0          1583
  NNC, novel donor and acceptor        1           0           1           180
  border junction                      0           1           0        154987



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
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_detection_nomulti.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn0
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn0_nomulti.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn2
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn2_nomulti.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn5
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn5_nomulti.csv"
     dbms=csv replace;
run;

proc export data=junc_count_apn10
     outfile="!MCLAB/event_Analysis/analysis_output/junctions_pacbio_vs_catalog_vs_star_support_apn10_nomulti.csv"
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



