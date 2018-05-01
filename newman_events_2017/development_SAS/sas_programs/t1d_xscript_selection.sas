ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splicing '/mnt/store/splice';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* Assign fusion detection to fragment detection so I can get to transcript lists quicker */

data fus_flags;
  set con.fusion_on_ge_apn5_v2;
  keep fusion_id flag_cd19_on flag_cd4_on flag_cd8_on;
run;

data frag2fus;
  set hg19.hg19_aceview_exon_fragment_info;
  keep fragment_id fusion_id;
run;

proc sort data=fus_flags;
  by fusion_id;
proc sort data=frag2fus;
  by fusion_id;
run;

data frag_flags;
  merge frag2fus (in=in1) fus_flags (in=in2);
  by fusion_id;
  if in1 and in2;
  drop fusion_id;
run;

data fragment_on_flags;
  set frag_flags;
  *if flag_cd19_on=. then flag_cd19_on=1;
  *if flag_cd4_on=. then flag_cd4_on=1;
  *if flag_cd8_on=. then flag_cd8_on=1;
   rename fragment_id=feature_id;
run;

data junc_on_flags;
   set event.t1d_flag_splicing_on_apn5;
   drop flag_splicing_all_ge5;
   rename event_id=feature_id;
run;

data feat_on_flags;
   set junc_on_flags fragment_on_flags;
run;

data feat2xs;
  set event.hg19_feature2xs2gene_exp_only2;
  drop gene_id;
run;

/* For each transcript, calculate the proportion of events detected (for each cell type) */

%macro calcProp(celltype);

data dtct_events;
   set feat_on_flags;
   where flag_&celltype._on=1;
   keep feature_id;
run;

proc sort data=feat2xs;
   by feature_id;
proc sort data=dtct_events;
   by feature_id;
run;

data dtct_feat2xs;
  merge feat2xs (in=in1) dtct_events (in=in2);
  by feature_id;
  if in2 then flag_feature_on=1; else flag_feature_on=0;
  flag_feature=1;
  if in1 then output;
run;

proc sort data=dtct_feat2xs;
   by transcript_id;
proc means data=dtct_feat2xs noprint;
   by transcript_id;
   var  flag_feature flag_feature_on;
   output out=num_feat_on_by_xs sum(flag_feature)=num_features sum(flag_feature_on)=num_features_dtct;
run;

data perc_feat_on_&celltype.;
  set num_feat_on_by_xs;
  perc_features_dtct=num_features_dtct/num_features;
  drop _TYPE_ _FREQ_;
run;

%mend;

%calcProp(cd4);
%calcProp(cd8);
%calcProp(cd19);

data perc_feat_on_cd4_2;
   set perc_feat_on_cd4;
   rename perc_features_dtct=perc_features_dtct_CD4
          num_features_dtct=num_features_dtct_CD4
	  num_features=num_features_CD4;
run;

data perc_feat_on_cd8_2;
   set perc_feat_on_cd8;
   rename perc_features_dtct=perc_features_dtct_CD8
          num_features_dtct=num_features_dtct_CD8
	  num_features=num_features_CD8;
run;

data perc_feat_on_cd19_2;
   set perc_feat_on_cd19;
   rename perc_features_dtct=perc_features_dtct_CD19
          num_features_dtct=num_features_dtct_CD19
	  num_features=num_features_CD19;
run;

proc sort data=perc_feat_on_cd4_2;
   by transcript_id;
proc sort data=perc_feat_on_cd8_2;
   by transcript_id;
proc sort data=perc_feat_on_cd19_2;
   by transcript_id;
run;

data perc_Feat_on_3cell;
  merge perc_feat_on_Cd4_2 perc_feat_on_Cd8_2 perc_feat_on_Cd19_2;
   by transcript_id;
run;

data flag_xs_by_perc_dtct;
   set perC_feat_on_3cell;
   if perc_features_dtct_CD4 >= 0.75 then flag_cd4_75perc=1; else flag_cd4_75perc=0;
   if perc_features_dtct_CD4 >= 0.5 then flag_cd4_50perc=1; else flag_cd4_50perc=0;
   if perc_features_dtct_CD8 >= 0.75 then flag_cd8_75perc=1; else flag_cd8_75perc=0;
   if perc_features_dtct_CD8 >= 0.5 then flag_cd8_50perc=1; else flag_cd8_50perc=0;
   if perc_features_dtct_CD19 >= 0.75 then flag_CD19_75perc=1; else flag_CD19_75perc=0;
   if perc_features_dtct_CD19 >= 0.5 then flag_CD19_50perc=1; else flag_CD19_50perc=0;
run;

proc freq noprint data=flag_xs_by_perc_dtct;
  tables flag_cd4_75perc*flag_cd8_75perc*flag_cd19_75perc / out=xs_count;
proc print data=xs_count;
run;
/*
                            flag_
 flag_cd4_    flag_cd8_     CD19_
   75perc       75perc     75perc    COUNT    PERCENT

     0            0           0      57445    66.5235
     0            0           1       2379     2.7550
     0            1           0        763     0.8836
     0            1           1        664     0.7689
     1            0           0        409     0.4736
     1            0           1        234     0.2710
     1            1           0       2112     2.4458
     1            1           1      22347    25.8787

total possible=22347+2112+234+409+664+763+2379 = 28908
*/

proc freq noprint data=flag_xs_by_perc_dtct;
  tables flag_cd4_50perc*flag_cd8_50perc*flag_cd19_50perc / out=xs_count;
proc print data=xs_count;
run;


/*
                            flag_
 flag_cd4_    flag_cd8_     CD19_
   50perc       50perc     50perc    COUNT

     0            0           0      47446
     0            0           1       2590
     0            1           0        809
     0            1           1        601
     1            0           0        489
     1            0           1        310
     1            1           0       2318
     1            1           1      31790

total possible=31790+2318+310+489+601+809+2590 = 38907

use this one for now...
*/


data xs2keep;
   set flag_xs_by_perc_dtct;
   where flag_cd4_50perc=1 or flag_cd8_50perc=1 or flag_cd19_50perc=1;
   keep transcript_id;
run;




/* Export transcript list */

data xs2av;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep gene_id transcript_id transcript_fasta_id;
run;

proc sort data=xs2keep;
  by transcript_id;
proc sort data=xs2av;
  by transcript_id;
run;

data xs2av_keep;
  merge xs2av (in=in1) xs2keep (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data xs_list;
  set xs2av_keep;
  keep transcript_fasta_id;
run;

data gene2xs_list;
   set xs2av_keep;
   keep gene_id transcript_fasta_id;
run;

proc export data=xs_list outfile="!MCLAB/event_analysis/references/t1d_xscripts_gene_exp2.txt"
   dbms=tab replace; putnames=no;
run;

proc export data=gene2xs_list outfile="!MCLAB/event_analysis/references/t1d_xscripts_gene_exp2_gene2xs.txt"
   dbms=tab replace; putnames=no;
run;


