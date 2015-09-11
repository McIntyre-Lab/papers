/* Junctions summary by gene - part 1: cat exons of events */

/* set libraries */
libname sgrloc '/media/jrbnewman/SAS_WRK1/sugrue/sas_data';
libname splice '/media/jrbnewman/SAS_WRK1/splice';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Cat splicing events on gene */
/* I think I might want to get the number of splicing events "detectable" (>10Reads, 1 sample) and "flag_all_on=1" and get these summaries

Table should look like:
Gene_id, cat_fusions, cat_exons, num_events_total, num_events_detectable, num_events_detect_all, num_events_de, num_intron_retention, num_IR_detectable, num_IR_detect_all, num_IR_de

*/

* get splicing events with data;

data splicing_de_summary;
  set sgrloc.results_by_jnct_fdr;
  if num_samples_exp ne 0;
  if feature1_id='intron' then feature1_id='';
  if feature2_id='intron' then feature2_id='';
  keep fusion_id  gene_id feature1_id feature2_id flag_intron_retention flag_all_on flag_p05 num_samples_exp;
run;

* Split file into (1) fusions and exons (2) flag_intron_retention num_samples_exp (3) flag_intron_retention flag_p05;

data events_for_cat_detect;
    set splicing_de_summary;
    keep gene_id fusion_id;
run;

data events_for_cat_de;
    set splicing_de_summary;
    if flag_p05=1 then output;
    keep gene_id fusion_id;
run;

proc sort data=events_for_cat_detect;
   by gene_id fusion_id;
run;

proc sort data=events_for_cat_de;
   by gene_id fusion_id;
run;


proc freq noprint data=events_for_cat_de;
   tables gene_id / out=events_cnt;
run;

proc sort data=events_cnt;
  by descending count;
run;
*max=30 DE events per gene;

proc sort data=events_for_cat_de;
   by gene_id fusion_id;
run;

data events_for_cat_de2; 
  array events[26] $ 17;

  retain events1-events26;

  set events_for_cat_de;
  by gene_id;
  
  if first.gene_id then do;
     call missing(of events1-events26);
     records = 0;
  end;

  records + 1;
  events[records]=fusion_id;
  if last.gene_id then output;
run;

  *clean up the output file;

data events_for_cat_de3;
  set events_for_cat_de2;
  length events_de_cat $ 468;
         events_de_cat= catx('|', OF events1-events30);
  rename records=num_events_de;
  drop events1-events26 fusion_id ;
  run;


* cat events detected;
proc freq noprint data=events_for_cat_detect;
   tables gene_id / out=events_cnt;
run;

proc sort data=events_cnt;
  by descending count;
run;
*max=93 detected events per gene;

proc sort data=events_for_cat_detect;
   by gene_id fusion_id;
run;

data events_for_cat_detect2; 
  array events[93] $ 17;

  retain events1-events93;

  set events_for_cat_detect;
  by gene_id;
  
  if first.gene_id then do;
     call missing(of events1-events93);
     records = 0;
  end;

  records + 1;
  events[records]=fusion_id;
  if last.gene_id then output;
run;

  *clean up the output file;

data events_for_cat_detect3;
  set events_for_cat_detect2;
  length events_detect_cat $ 1674;
         events_detect_cat= catx('|', OF events1-events93);
  rename records=num_events_detected;
  drop events1-events93 fusion_id ;
  run;

* merge together;
proc sort data=events_for_cat_detect3;
   by gene_id;
run;

proc sort data=events_for_cat_de3;
   by gene_id;
run;

proc sort data=exons_for_cat3;
   by gene_id;
run;

data exons_events_de oops1 oops2;
   merge events_for_cat_de3 (in=in1) exons_for_cat3 (in=in2);
   by gene_id;
   if in1 and in2 then output exons_events_de; 
   else if in1 then output oops1;
   else output oops2;
run;


data exons_events_de2;
   merge exons_events_de events_for_cat_detect3;
   by gene_id;
   if events_de_cat='' then events_de_cat='.';
   if exons_se_cat='' then exons_se_cat='.';
   if events_detect_cat='' then events_detect_cat='.';
run;


/* Make permenant */

data sugrue.exons_splicing_events_de;
    set exons_events_de2;
run;




