ods listing; ods html;

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/concannon/useful_human_data/aceview_hg19/fusions/sas_data';

/* Need to check:
T1D vs AI: number of possible events detected
*/

/* Merge cleaned event data with gene lists */

data events_w_flags;
   set splicing.flag_splicing_by_gene_dtct_v2;
   drop gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on;
run;

data all_events;
  set splice.splicing_Events_annotations;
  keep event_id gene_id;
run;


proc sort data=events_w_flags;
   by event_id ;
proc sort data=all_events;
   by event_id ;
run;

data all_events_w_flags;
  merge all_events (in=in1) events_w_flags (in=in2);
  by event_id;
  if not in1 then do;
     flag_cd4_on=0; flag_cd8_on=0; flag_cd19_on=0;
     end;
run;


data gene_flags;
  set con.flag_gene_detection_by_cell;
  drop num_exons_exp_: ;
run;

proc sort data=all_events_w_flags;
  by gene_id;
proc sort data=gene_flags; 
  by gene_id;
run;

data all_events_w_geneflags;
  merge all_events_w_flags (in=in1) gene_flags (in=in2);
  by gene_id;
  if in1 and in2;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d nodup;
   by gene_id;
proc sort data=all_events_w_geneflags;
   by gene_id;
run;

data events_w_flags_immuno;
   merge all_events_w_geneflags (in=in1) ai (in=in2) t1d (in=in3);
   by gene_id;
   if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;


data flag_event_any_on;
  set events_w_flags_immuno;
  if sum(flag_cd4_gene_on,flag_cd8_gene_on,flag_cd19_gene_on) > 1 then do;
      if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) > 0 then flag_event_cell_any_on=1;
      else flag_event_cell_any_on=0;
      end;
  else flag_event_cell_any_on=. ;
run;

proc freq data=flag_event_any_on;
  where flag_immuno_gene=1;
  tables flag_event_cell_any_on*flag_immunobase_diabetes_gene;
run;


/*
 flag_event_cell_any_on
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 496696 | 142507 | 639203
          |  74.73 |  21.44 |  96.17
          |  77.71 |  22.29 |
          |  96.20 |  96.05 |
 ---------+--------+--------+
        1 |  19609 |   5861 |  25470
          |   2.95 |   0.88 |   3.83
          |  76.99 |  23.01 |
          |   3.80 |   3.95 |
 ---------+--------+--------+
 Total      516305   148368   664673
             77.68    22.32   100.00
*/

/* Count genes */

proc means data=flag_event_any_on noprint;
   by gene_id;
   var flag_event_cell_any_on flag_immuno_gene flag_immunobase_diabetes_gene;
   output out=genes_any_event_on max=;
run;

proc freq data=genes_any_event_on;
   where flag_immuno_gene=1;
   tables flag_event_cell_any_on*flag_immunobase_diabetes_gene;
run;

/*
   flag_event_cell_any_on
             flag_immunobase_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |    114 |     46 |    160
            |   6.97 |   2.81 |   9.78
            |  71.25 |  28.75 |
            |   9.18 |  11.68 |
   ---------+--------+--------+
          1 |   1128 |    348 |   1476
            |  68.95 |  21.27 |  90.22
            |  76.42 |  23.58 |
            |  90.82 |  88.32 |
   ---------+--------+--------+
   Total        1242      394     1636
               75.92    24.08   100.00
*/
