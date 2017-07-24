ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Divide BLAST hits into annotated and unannotated events (features). For annotated features,
   I want to plot the APN distribution (density plot) of annotated features with a hit
     vs annotated features without a hit */

data features_w_hit;
   set event.pacbio_hit_correct_gene;
   keep feature_id;
run;

data annot_feat;
   set event.features_w_annotations_nomulti;
   where flag_feature_on=1;
   if feature_type="fragment" then output;
   if feature_type="splicing" and flag_junction_annotated=1 then output;
   keep feature_id feature_type flag_singleton;
run;

proc sort data=features_w_hit nodup;
  by feature_id;
proc sort data=annot_feat;
  by feature_id;
run;

data flag_annot_hit;
  merge annot_feat (in=in1) features_w_hit (in=in2);
  by feature_id;
  if in2 then flag_pacbio_hit=1; else flag_pacbio_hit=0;
run;

proc freq data=flag_annot_hit;
tables flag_pacbio_hit;
run;

proc freq data=flag_annot_hit noprint;
 tables flag_pacbio_hit*flag_singleton*feature_type / out=hits_by_type;
proc print data=hits_by_type;
run;

/*
                                                 Cumulative    Cumulative
     flag_pacbio_hit    Frequency     Percent     Frequency      Percent
     --------------------------------------------------------------------
                   0      174133       72.62        174133        72.62
                   1       65659       27.38        239792       100.00

   flag_
  pacbio_      flag_      feature_
    hit      singleton      type      COUNT    PERCENT

     0           0        fragment    40264    16.7912
     0           0        splicing    47337    19.7409
     0           1        fragment    86532    36.0863
     1           0        fragment     3796     1.5830
     1           0        splicing    42475    17.7133
     1           1        fragment    19388     8.0853



*/


/* Now calc by gene */

data feat2gene;
   set event.feature2xs2gene_exp_only_nomulti;
   keep gene_id feature_id;
run;

proc sort data=feat2gene nodup;
   by feature_id gene_id;
proc sort data=flaG_annot_hit;
   by feature_id;
run;

data annot_hit_w_gene oops;
   merge flag_annot_hit (in=in1) feat2gene (in=in2);
   by feature_id;
   if in1 and in2 then output annot_hit_w_gene;
   else if in1 then output oops;
run;


data gene_hit gene_nohit;
   set annot_hit_w_gene;
   if flag_pacbio_hit=1 then output gene_hit;
   else output gene_nohit;
   keep gene_id;
run;

proc sort data=gene_hit nodup;
   by gene_id;
proc sort data=gene_nohit nodup;
   by gene_id;
run;

data genes_w_hit;
  merge gene_hit (in=in1) gene_nohit (in=in2);
  by gene_id;
  if in1 then flag_pacbio_hit=1; else flag_pacbio_hit=0;
run;

proc freq data=genes_w_hit;
  tables flag_pacbio_hit;
run;

/*
                                                Cumulative    Cumulative
    flag_pacbio_hit    Frequency     Percent     Frequency      Percent
    --------------------------------------------------------------------
                  0       16461       78.86         16461        78.86
                  1        4414       21.14         20875       100.00
*/


/* Calc mean APN per feature, merge into list here */

data frag_counts;
    set event.mm10_refseq_fragment_counts;
    where sample_id ? "NSC";
    keep sample_id fragment_id apn;
    rename fragment_id=feature_id;
run;

data event_counts;
    set event.mm10_refseq_splicing_counts;
    where sample_id ? "NSC";
    keep sample_id event_id apn;
    rename event_id=feature_id;
run;

data feature_counts;
   set event_counts frag_counts;
run;

proc sort data=feature_counts;
   by feature_id sample_id;
proc means data=feature_counts noprint;
   by feature_id;
   var apn;
   output out=mean_apn_by_feat mean=mean_apn_npc;
run;

proc sort data=mean_apn_by_feat;
  by feature_id;
proc sort data=flag_annot_hit;
   by feature_id;
run;

data annot_hit_w_apn;
  merge flag_annot_hit (in=in1) mean_apn_by_feat (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* before I export, I want to check the distribution of APN here */

proc sort data=annot_hit_w_apn;
  by flag_pacbio_hit feature_type;
proc means data=annot_hit_w_apn noprint;
  by flag_pacbio_hit feature_type;
  var mean_apn_npc;
  output out=apn_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc print data=apn_distrib ; run;

/*
         flag_
        pacbio_    feature_
 Obs      hit        type      _TYPE_    _FREQ_            mean              sd             min

  1        0       fragment       0      126796    33.186075398    263.14217635    0.0001643115
  2        0       splicing       0       47337    15.180574393    122.30641567     0.024691358
  3        1       fragment       0       23184    104.35429151    353.93787073    0.0138888889
  4        1       splicing       0       42475    56.860959603    280.22984938    0.1728395062



 Obs              q1          median              q3             max

  1     0.3414634146     1.762245102    10.200922266    25968.547727
  2     0.3456790123     1.024691358    2.9382716049     4569.345679
  3     11.453003381    31.142351915    85.490197568     13274.41954
  4     4.3086419753    12.586419753    37.777777778    19137.549383

Okay, looks like the APN is generally higher for features with PB hits
*/

/* Make permenant and export for plots */

data event.pacbio_hits_annot_feat_mean_apn;
   set annot_hit_w_apn;
   drop _TYPE_ _FREQ_;
run;


proc export data=event.pacbio_hits_annot_feat_mean_apn
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_hits_annotated_features_apn.csv"
     dbms=csv replace;
run;


