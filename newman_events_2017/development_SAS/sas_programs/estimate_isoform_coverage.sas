
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname splicing '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Make a "master" chunk/junction to transcript database, flagging features that are unique to a single isoform */

data splicing_info;
   set splicing.splicing_events_annotations;
   where transcript_id ne ' ';
   keep event_id gene_id chr transcript_id num_transcripts;
   rename event_id=feature_id;
run;

data chunks_info;
   set conesa.exon_chunk_w_info;
   keep exonchunk_id chrom gene_id transcript_id transcript_count;
   rename transcript_count=num_transcripts chrom=chr;
   rename exonchunk_id=feature_id;
run;

data feature_info;
   length feature_id $2550.;
   length transcript_id $450.;
   set splicing_info chunks_info;
run;

* unstack transcript_id;

data feature_info_unstack;
   set feature_info;
   length transcript_id2 $35.; 
   do i=1 by 1 while(scan(transcript_id,i,'|') ^=' '); 
transcript_id2=scan(transcript_id,i,'|'); 
drop transcript_id i;
rename transcript_id2=transcript_id;
output; 
end; 
run;


/* Merge with counts */

data chunk_counts;
   set conesa.coverage_exon_chunks;
   where sample_id ? "NSC";
   keep sample_id apn fusion_id;
   rename fusion_id=feature_id;
run;

data splicing_counts;
   set conesa.coverage_splicing;
   where sample_id ? "NSC";
   keep sample_id apn fusion_id;
   rename fusion_id=feature_id;
run;

data all_counts;
   length feature_id $2550.;
   set  splicing_counts chunk_counts;
run;

/* Flag transcripts if they have at least one unique piece */

data uniq_xs;
   set feature_info_unstack;
   if num_transcripts=1;
   keep transcript_id;
run;

data all_xs;
   set feature_info_unstack;
   keep transcript_id;
run;

proc sort data=uniq_xs nodup;
   by transcript_id;
proc sort data=all_xs nodup;
   by transcript_id;
run;

* 14675 of 16104 isoforms with at least one unique piece ;


data flag_xs_w_uniq;
   merge all_xs (in=in1) uniq_xs (in=in2);
   by transcript_id;
   if in2 then flag_xscript_uniq_feature=1;
   else flag_xscript_uniq_feature=0;
run;

/* Flag features if on/off */

data flag_on_counts;
   set all_counts;
   if apn > 0 then flag_feature_on=1;
   else flag_feature_on=0;
run;

proc sort data=flag_on_counts;
   by feature_id;
proc means data=flag_on_counts noprint;
   by feature_id;
   var flag_feature_on;
   output out=flag_feature_on mean=perc_on;
run;

data flag_feature_on2;
   set flag_feature_on;
   if perc_on lt 0.5 then flag_feature_on=0;
   else flag_feature_on=1;
run;

* merge in feature2xs;

proc sort data=flag_feature_on2;
   by feature_id;
proc sort data=feature_info_unstack;
   by feature_id;
run;

data feature_info_w_flag oops1;
  merge feature_info_unstack (in=in1) flag_feature_on2 (in=in2);
  by feature_id;
  if in1 and in2 then output feature_info_w_flag;
  else if in2 then output oops1;
run;

/* Count transcripts with at least one unique piece that is detected */

data xs_w_on_uniq;
   set feature_info_w_flag;
   if num_transcripts=1 and flag_feature_on=1;
   keep transcript_id;
run;

proc sort data=xs_w_on_uniq nodup;
   by transcript_id;
run;

* 13159 of 16104 isoforms with at least one unique piece detected ;

/* For transcripts with at least one unique piece (on or otherwise), I want to average across uniq pieces and then all pieces */

/* average counts */

proc sort data=all_counts;
   by feature_id;
proc means data=all_counts noprint;
  by feature_id;
  var apn;
  output out=mean_counts mean=;
run;


proc sort data=feature_info_w_flag;
   by feature_id;
proc sort data=mean_counts;
   by feature_id;
run;

data xs_w_counts;
  merge mean_counts (in=in1) feature_info_w_flag (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* Subset only unique pieces and estimate abundance -- unique counts per isoform */

data uniq;
 set xs_w_counts;
 where num_transcripts=1;
run;

proc sort data=uniq;
   by transcript_id;
proc means data=uniq noprint;
   by transcript_id;
   var apn;
   output out=estimate_uniq mean=;
run;
*14675 xs;

/* Subset only unique "on" pieces and estimate abundance -- unique counts per isoform: only detected features */

data uniq_dtct;
 set xs_w_counts;
 where num_transcripts=1 and flag_feature_on=1;
run;

proc sort data=uniq_dtct;
   by transcript_id;
proc means data=uniq_dtct noprint;
   by transcript_id;
   var apn;
   output out=estimate_uniq_dtct mean=;
run;
*13159 xs;

/* Subset all pieces and estimate abundance -- all isoforms */

data all;
 set xs_w_counts;
run;


proc sort data=all;
   by transcript_id;
proc means data=all noprint;
   by transcript_id;
   var apn;
   output out=estimate_all mean=;
run;
*16104 xs;

/* Subset all pieces "on" pieces and estimate abundance -- all isoforms: only detected features */

data all_dtct;
 set xs_w_counts;
 where flag_feature_on=1;
run;

proc sort data=all_dtct;
   by transcript_id;
proc means data=all_dtct noprint;
   by transcript_id;
   var apn;
   output out=estimate_all_dtct mean=;
run;
*16066 xs;

/* Merge all together */

data estimate_uniq2;
  set estimate_uniq;
  drop _TYPE_ _FREQ_;
  rename apn=apn_unique;
run;

data estimate_uniq_dtct2;
  set estimate_uniq_dtct;
  drop _TYPE_ _FREQ_;
  rename apn=apn_unique_dtct;
run;

data estimate_all2;
  set estimate_all;
  drop _TYPE_ _FREQ_;
  rename apn=apn_all;
run;

data estimate_all_dtct2;
  set estimate_all_dtct;
  drop _TYPE_ _FREQ_;
  rename apn=apn_all_dtct;
run;

proc sort data=estimate_uniq2;
   by transcript_id;
proc sort data=estimate_uniq_dtct2;
   by transcript_id;
proc sort data=estimate_all2;
   by transcript_id;
proc sort data=estimate_all_dtct2;
   by transcript_id;
run;

data all_estimates;
   merge estimate_uniq2 estimate_uniq_dtct2 estimate_all2 estimate_all_dtct2;
   by transcript_id;
run;


* check distribution ;

proc means data=all_estimates noprint;
   var apn_unique;
   output out=apn_unq_distrib p5=p5 p95=p95 median=median min=min max=max q3=q3 q1=q1;
run;

/* Min=0
   P5=0
   Q1=0.7166666667
   MED=3.3055555556
   Q3=11.164170316
   P95=70.876376909
   MAX=35525.061728
*/

proc means data=all_estimates noprint;
   var apn_unique_dtct;
   output out=apn_unq_dtct_distrib p5=p5 p95=p95 median=median min=min max=max q3=q3 q1=q1;
run;

/* Min=0.001641497
   P5=0.3161764706
   Q1=1.5409713743
   MED=4.5669191919
   Q3=13.66607015
   P95=82.037777982
   MAX=35525.061728
*/

proc means data=all_estimates noprint;
   var apn_all;
   output out=apn_all_distrib p5=p5 p95=p95 median=median min=min max=max q3=q3 q1=q1;
run;


/* Min=0
   P5=1.4908171921
   Q1=8.057025045
   MED=22.482478057
   Q3=61.100259395
   P95=316.64016858
   MAX=10071.640032
*/
proc means data=all_estimates noprint;
   var apn_all_dtct;
   output out=apn_all_dtct_distrib p5=p5 p95=p95 median=median min=min max=max q3=q3 q1=q1;
run;


/* Min=0.0105162524
   P5=1.7523074191
   Q1=8.4894506369
   MED=23.695784513
   Q3=64.632000684
   P95=337.45869084
   MAX=10071.640032
*/



/* Flag if APN>5 */

data all_estimates_gt5;
  set all_estimates;
  if apn_unique ge 5 then flag_apn_unique_gt5=1; else flag_apn_unique_gt5=0;
  if apn_unique_dtct ge 5 then flag_apn_unique_dtct_gt5=1; else flag_apn_unique_dtct_gt5=0;
  if apn_all ge 5 then flag_apn_all_gt5=1; else flag_apn_all_gt5=0;
  if apn_all_dtct ge 5 then flag_apn_all_dtct_gt5=1; else flag_apn_all_dtct_gt5=0;
run;

*quick count of transcripts APN ge 5;

data uniq uniq_dtct all all_dtct;
   set all_estimates_gt5;
   if flag_apn_unique_gt5=1 then output uniq;
   if flag_apn_unique_dtct_gt5=1 then output uniq_dtct;
   if flag_apn_all_gt5=1 then output all;
   if flag_apn_all_dtct_gt5=1 then output all_dtct;
run;

/* COUNTS:
   Unique >5 APN		6019
   Unique detected >5 APN	6264
   All >5 APN			13464
   All detected >5 APN		13630

Plot the distribution of each in Python!!
*/


data check;
   set conesa.isoform_estimates;
   where apn_all_dtct ne .;
run;

