ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Match computed introns to intron retention events. Flag IR events that have no associated intron
   (i.e. one or both associated exons are associated with multigene exonic regions) */

data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;


data intron_coord;
  set mm10.mm10_introns_from_fusions_si;
  keep intron_id chr intron_start intron_stop gene_id;
run;
*227967 computed introns;


data ir_events;
   set evspl.splicing_events_annot_refseq;
   where flag_intron_retention=1;
   keep event_id gene_id chr strand feature1_start  feature1_stop
        feature2_start  feature2_stop;
run;
*247983 IR events;

proc sort data=intron_coord;
  by gene_id;
proc sort data=ir_events;
  by gene_id;
proc sort data=exp_genes;
  by gene_id;
run;

data intron_coord2;
  merge intron_coord (in=in1) exp_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

data ir_events2;
  merge ir_events (in=in1) exp_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;


proc sort data=intron_coord2;
   by chr intron_start intron_stop;
proc sort data=ir_events2;
   by chr feature2_start feature2_stop;
run;

/*
For + strand: merge intron_start on feature2_start (+1) 
For - strand: merge intron_stop on feature1_stop (+1)

Remember: intron coordinates in SAS dataset are for BED files (which are 0-based), so I need to add 1bp to shift it
          so that it now lines up with the IR event coordinates */

data pos_ir neg_ir;
   set ir_events2;
   if strand="+" then output pos_ir;
   else output neg_ir;
run;

data intron_coord3;
   set intron_coord2;
   feature1_stop=intron_stop+1;
   feature2_start=intron_start+1;
run;

proc sort data=pos_ir;
   by chr feature2_start;
proc sort data=intron_coord3;
   by chr feature2_start;
run;

data pos_ir_w_intron neg_intron;
  merge intron_coord3 (in=in1) pos_ir (in=in2);
  by chr feature2_start;
  if in1 and in2 then output pos_ir_w_intron;
  else if in1 then output neg_intron;
run;

proc sort data=neg_intron;
   by chr feature1_stop;
proc sort data=neg_ir;
   by chr feature1_stop;
run;

data neg_ir_w_intron no_ir;
  merge neg_intron (in=in1) neg_ir (in=in2);
  by chr feature1_stop;
  if in1 and in2 then output neg_ir_w_intron;
  else if in1 then output no_ir;
run;

data ir_all;
   set ir_events;
   keep event_id;
run;

data stack_ir_w_intron;
  set pos_ir_w_intron neg_ir_w_intron;
  keep event_id intron_id gene_id;
run;

proc sort data=stack_ir_w_intron;
   by event_id;
proc sort data=ir_all;
   by event_id;
run;

data flag_no_intron;
  merge ir_all (in=in1) stack_ir_w_intron (in=in2);
  by event_id;
  if in1 and in2 then flag_no_computed_intron=0;
  else if in1 then flag_no_computed_intron=1;
  else flag_no_computed_intron=.; *check to make sure we haven't added any extra IR somehow;
run;


proc freq data=flag_no_intron;
  tables flag_no_computed_intron;
run;

/*

         flag_no_                             Cumulative    Cumulative
  computed_intron    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      154934       62.48        154934        62.48
                1       93049       37.52        247983       100.00

*/

/* Make permenant */

data event.ir_events_w_intron;
   set flag_no_intron;
run;

