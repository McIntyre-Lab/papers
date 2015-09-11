/* Junctions summary by gene - part 2: count DE splicing events */

/* set libraries */
libname sgrloc '/mnt/SAS_WRK1/sugrue/sas_data';
libname splice '/mnt/data/splice';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Count splicing events and intron_retention events detected */

* get splicing events with data;

data splicing_de_summary;
  set sgrloc.results_by_jnct_fdr;
  if num_samples_exp ne 0;
  if feature1_id='intron' then feature1_id='';
  if feature2_id='intron' then feature2_id='';
  keep fusion_id gene_id feature1_id feature2_id flag_intron_retention flag_all_on flag_p05 num_samples_exp;
run;

proc sort data=splicing_de_summary;
   by gene_id;
run;

data detected_events;
    set splicing_de_summary;
    event_count=1;
    keep event_count event_id gene_id flag_intron_retention;
run;

proc means data=detected_events noprint;
    by gene_id;
    var event_count flag_intron_retention;
output out=detected_events_sum sum=;
run;


/* Count splicing events and intron_retention events expressed in both groups */


data analyzed_events;
    set splicing_de_summary;
    if flag_all_on=1;
    event_count=1;
    keep event_count fusion_id gene_id flag_intron_retention;
run;

proc means data=analyzed_events noprint;
    by gene_id;
    var event_count flag_intron_retention;
output out=analyzed_events_sum sum=;
run;

/* Count splicing events and intron_retention events differentially expressed */


data diffexp_events;
    set splicing_de_summary;
    if flag_p05=1;
    event_count=1;
    keep event_count fusion_id gene_id flag_intron_retention;
run;

proc means data=diffexp_events noprint;
    by gene_id;
    var event_count flag_intron_retention;
output out=diffexp_events_sum sum=;
run;

/* Merge datasets into one big table */


* se=splicing event, ir=intron retention;

data detected_events_sum_2;
   set detected_events_sum;
   rename 
          event_count=detected_se_cnt
          flag_intron_retention=detected_ir_cnt;
run;


data analyzed_events_sum_2;
   set analyzed_events_sum;
   rename 
          event_count=analyzed_se_cnt
          flag_intron_retention=analyzed_ir_cnt;
run;

data diffexp_events_sum_2;
   set diffexp_events_sum;
   rename 
          event_count=diffexp_se_cnt
          flag_intron_retention=diffexp_ir_cnt;
run;

* sort datasets;

proc sort data=detected_events_sum_2;
   by gene_id;
run;


proc sort data=analyzed_events_sum_2;
   by gene_id;
run;


proc sort data=diffexp_events_sum_2;
   by gene_id;
run;

* merge datasets;

* detected and analyzed splicing events;
data se_summary_1 oops;
   merge detected_events_sum_2 (in=in1) analyzed_events_sum_2 (in=in2);
   by gene_id;
   if in1 and in2 then output se_summary_1;
   else if in1 then do;
       analyzed_se_cnt=0;
       analyzed_ir_cnt=0;
       output se_summary_1;
       end;
   else if in2 then output oops;
run;


proc sort data=se_summary_1;
    by gene_id;
run;

* detected, analyzed and DE splicing events;
data se_summary_2 oops;
   merge se_summary_1 (in=in1) diffexp_events_sum_2 (in=in2);
   by gene_id;
   if in1 and in2 then output se_summary_2;
   else if in1 then do;
       diffexp_se_cnt=0;
       diffexp_ir_cnt=0;
        output se_summary_2;
       end;
   else if in2 then output oops;
run;

proc sort data=se_summary_2;
    by gene_id;
run;

/* Make permenant */

data sugrue.splicing_by_gene_summary;
set se_summary_2;
run;
