ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Merge on-flags, mean APN, introns, fragments and fusions.
   Then calculate the ratio between IR event and nearest fusion and fragment and export for plots

   I want to plot density plots of IR to fusion/fragment ratios and scatter plots of
   IR apn to fusion/fragment APN

   If IR event not detected, then we should drop the intron check */

/* Get detection flags */
data int_on;
    set event.flag_intron_on;
    keep intron_id flag_intron_nsc_on;
run;

data ir_on;
   set event.splicing_on_apn_gt0;
   keep event_id flag_splicing_on;
run; 

data fus_on;
   set event.flag_fusion_on;
   keep fusion_id flag_fusion_nsc_on;
run;

data frag_on;
   set event.fragments_on_apn_gt0;
   keep fragment_id flag_fragment_on;
run;

/* Get mean APN */
data mean_int_apn;
   set event.mean_apn_intron_nsc;
run;

data mean_frag_apn;
   set event.mean_apn_fragment_nsc;
run;

data mean_fus_apn;
   set event.mean_apn_fusion_nsc;
run;

data mean_IR_apn;
   set event.mean_apn_IR_nsc;
run;

/* Get intron-to-IR and intron-to-fusion/fragment and merge */
data intron_to_event;
  set event.ir_events_w_intron;
   where flag_no_computed_intron=0; *only keep those with computed introns;
run;

data intron_to_fus_frag;
   set mm10.mm10_introns_to_fragment_fusion;
   drop chr intron_start intron_stop;
run;

proc sort data=intron_to_event;
   by intron_id;
proc sort data=intron_to_fus_frag;
   by intron_id;
run;

data intron_to_ir_fus_frag;
  merge intron_to_event (in=in1) intron_to_fus_frag (in=in2);
  by intron_id;
  if in1 and in2;
run;

/* Merge all together */
proc sort data=int_on;
  by intron_id;
proc sort data=mean_int_apn;
  by intron_id;
run;
data mean_int_apn_w_dtct;
   merge int_on (in=in1) mean_int_apn (in=in2);
   by intron_id;
   if in1 and in2;
run;

proc sort data=ir_on;
  by event_id;
proc sort data=mean_ir_apn;
  by event_id;
run;
data mean_ir_apn_w_dtct;
   merge ir_on (in=in1) mean_ir_apn (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=fus_on;
  by fusion_id;
proc sort data=mean_fus_apn;
  by fusion_id;
run;
data mean_fus_apn_w_dtct;
   merge fus_on (in=in1) mean_fus_apn (in=in2);
   by fusion_id;
   if in1 and in2;
run;

proc sort data=frag_on;
  by fragment_id;
proc sort data=mean_frag_apn;
  by fragment_id;
run;
data mean_frag_apn_w_dtct;
   merge frag_on (in=in1) mean_frag_apn (in=in2);
   by fragment_id;
   if in1 and in2;
run;

data frag_5prime_mean_dtct;
   set mean_frag_apn_w_dtct;
   rename fragment_id=fragment_id_5prime flag_fragment_on=flag_fragment_on_5prime
          mean_apn_fragment=mean_apn_fragment_5prime;
run;

data frag_3prime_mean_dtct;
   set mean_frag_apn_w_dtct;
   rename fragment_id=fragment_id_3prime flag_fragment_on=flag_fragment_on_3prime
          mean_apn_fragment=mean_apn_fragment_3prime;
run;

data fus_5prime_mean_dtct;
   set mean_fus_apn_w_dtct;
   rename fusion_id=fusion_id_5prime flag_fusion_nsc_on=flag_fusion_on_5prime
          mean_apn_fusion=mean_apn_fusion_5prime;
run;

data fus_3prime_mean_dtct;
   set mean_fus_apn_w_dtct;
   rename fusion_id=fusion_id_3prime flag_fusion_nsc_on=flag_fusion_on_3prime
          mean_apn_fusion=mean_apn_fusion_3prime;
run;
   
proc sort data=frag_5prime_mean_dtct;
   by fragment_id_5prime;
proc sort data=frag_3prime_mean_dtct;
   by fragment_id_3prime;
proc sort data=fus_5prime_mean_dtct;
   by fusion_id_5prime;
proc sort data=fus_3prime_mean_dtct;
   by fusion_id_3prime;
proc sort data=mean_ir_apn_w_dtct;
   by event_id;
proc sort data=mean_int_apn_w_dtct;
   by intron_id;
proc sort data=intron_to_ir_fus_frag;
   by intron_id;
run;

data intron2ir_p1;
  merge intron_to_ir_fus_frag (in=in1) mean_int_apn_w_dtct (in=in2);
  by intron_id;
  if in1 and in2;
run;

proc sort data=intron2ir_p1;
   by event_id;
run;

data intron2ir_p2;
   merge intron2ir_p1 (in=in1) mean_ir_apn_w_dtct (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=intron2ir_p2;
  by fusion_id_5prime;
run;

data intron2ir_p3;
  merge intron2ir_p2 (in=in1) fus_5prime_mean_dtct (in=in2);
  by fusion_id_5prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p3;
  by fusion_id_3prime;
run;

data intron2ir_p4;
  merge intron2ir_p3 (in=in1) fus_3prime_mean_dtct (in=in2);
  by fusion_id_3prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p4;
  by fragment_id_5prime;
run;

data intron2ir_p5;
  merge intron2ir_p4 (in=in1) frag_5prime_mean_dtct (in=in2);
  by fragment_id_5prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p5;
  by fragment_id_3prime;
run;

data intron2ir_p6;
  merge intron2ir_p5 (in=in1) frag_3prime_mean_dtct (in=in2);
  by fragment_id_3prime;
  if in1 and in2;
run;

/* Calc IR ratios */
data ir_ratios;
  set intron2ir_p6;
  ir_to_5fus_ratio=mean_apn_ir/mean_apn_fusion_5prime;
  ir_to_5frag_ratio=mean_apn_ir/mean_apn_fragment_5prime;
  ir_to_3fus_ratio=mean_apn_ir/mean_apn_fusion_3prime;
  ir_to_3frag_ratio=mean_apn_ir/mean_apn_fragment_3prime;
  mean_apn_fusion_mean=(mean_apn_fusion_5prime+mean_apn_fusion_3prime)/2;
  mean_apn_fragment_mean=(mean_apn_fragment_5prime+mean_apn_fragment_3prime)/2;
  ir_to_mean_fus_ratio=mean_apn_ir/mean_apn_fusion_mean;
  ir_to_mean_frag_ratio=mean_apn_ir/mean_apn_fragment_mean;
run;

/* Make permenant */

data event.intron2ir_frag_fus_all_info;
   set ir_ratios;
run;

/* Flag if in expressed gene and export */

data exp_gene;
  set event.flag_gene_expressed;
  keep gene_id flag_gene_expressed;
run;

proc sort data=exp_gene;
  by gene_id;
proc sort data=ir_ratios;
  by gene_id;
run;

data ir_ratios_w_gene;
  merge ir_ratios (in=in1) exp_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc export data=ir_ratios_w_gene outfile="!MCLAB/event_analysis/analysis_output/event_analysis_intron_ir_analysis_nomulti.csv"
   dbms=csv replace;
run;

