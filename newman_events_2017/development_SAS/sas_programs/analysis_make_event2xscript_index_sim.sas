/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* I am treating each test (test1 and test2) as separate "conditions", and then each simulation (1,2,3) as replicates
   for each simulation. For each test then, I need to identify the set of genes that are expressed, skipping
   genes with multigene regions (fusions)

   First, for each test, I need to flag fusions as on/off (APN>0) then use this information
   to create a filtered event2transcript2gene index */


data fus_counts;
   set event.mm10_refseq_fusion_counts_sim;
   length test $5.;
   test=scan(sample_id,2,"_");
   if apn > 0 then flag_fusion_gt0=1; else flag_fusion_gt0=0;
run;

proc sort data=fus_counts;
   by fusion_id test;
proc means data=fus_counts noprint;
   by fusion_id test;
   var flag_fusion_gt0;
   output out=perc_trt_on mean=;
run;

data test1 test2;
   set perc_trt_on;
   if test="test1" then output test1;
   else if test="test2" then output test2;
run;

data test1_2;
  set test1;
  if flag_fusion_gt0 ge 0.5 then flag_fusion_test1_apn0=1;
  else flag_fusion_test1_apn0=0;
  keep fusion_id flag_fusion_test1_apn0;
run;

data test2_2;
  set test2;
  if flag_fusion_gt0 ge 0.5 then flag_fusion_test2_apn0=1;
  else flag_fusion_test2_apn0=0;
  keep fusion_id flag_fusion_test2_apn0;
run;

proc sort data=test1_2;
  by fusion_id;
proc sort data=test2_2;
  by fusion_id;
run;

data event.flag_fusion_on_apn0_sim;
   merge test1_2 (in=in1) test2_2 (in=in2);
   by fusion_id;
   if in1 and in2;
run;

/* check overlap */

proc freq data=event.flag_fusion_on_apn0_sim;
   tables flag_fusion_test1_apn0*flag_fusion_test2_apn0;
run;

/*

   flag_fusion_test1_apn0
             flag_fusion_test2_apn0

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 | 101556 |  30764 | 132320
            |  36.91 |  11.18 |  48.09
            |  76.75 |  23.25 |
            |  77.67 |  21.30 |
   ---------+--------+--------+
          1 |  29205 | 113652 | 142857
            |  10.61 |  41.30 |  51.91
            |  20.44 |  79.56 |
            |  22.33 |  78.70 |
   ---------+--------+--------+
   Total      130761   144416   275177
               47.52    52.48   100.00
*/

/* Now make simulation-specific event2xs indices */

* Get genes that do not have multigene fusions;

data genes_w_mult;
  set mm10.mm10_refseq_fusion_si_info_v2;
  if flag_multigene=1 then output genes_w_mult;
  keep primary_gene_id;
  rename primary_gene_id=gene_id;
run;

data event2xs2gene;
  set event.feature2xs2gene;
run;

proc sort data=genes_w_mult nodup;
  by gene_id; *10066 genes to remove;
proc sort data=event2xs2gene;
  by gene_id;
run;

data event2xs2gene_nomulti;
  merge event2xs2gene (in=in1) genes_w_mult (in=in2);
  by gene_id;
  if in2 then delete;
  else output;
run;

data fus2gene;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep primary_gene_id fusion_id;
  rename primary_gene_id=gene_id;
run;

proc sort data=fus2gene nodup;
  by fusion_id gene_id;
run;

%macro eventIndex(test);

data fus_on;
  set event.flag_fusion_on_apn0_sim;
  keep fusion_id flag_fusion_&test._apn0;
run;

proc sort data=fus_on;
  by fusion_id;
run;

data fus2gene_on;
  merge fus_on (in=in1) fus2gene (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2gene_on;
   by gene_id;
proc means data=fus2gene_on noprint;
   by gene_id;
   var flag_fusion_&test._apn0;
   output out=fus_on_per_gene sum=fusions_dtct;
run;

data genes_exp;
   set fus_on_per_gene;
   where fusions_dtct > 0;
   keep gene_id;
run;

proc sort data=genes_exp;
  by gene_id;
proc sort data=event2xs2gene_nomulti;
  by gene_id;
run;

data event2xs2gene_exp_only_nomulti;
  merge event2xs2gene_nomulti (in=in1) genes_exp (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Make permenant */
data event.feature2xs2gene_exp_nomult_&test.;
   set event2xs2gene_exp_only_nomulti;
run;

%mend;

%eventIndex(test1);
%eventIndex(test2);


