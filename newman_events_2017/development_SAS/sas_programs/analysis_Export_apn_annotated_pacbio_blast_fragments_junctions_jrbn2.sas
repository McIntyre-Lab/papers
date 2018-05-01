ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Divide BLAST hits into annotated and unannotated events (features). For annotated features,
   I want to plot the APN distribution (density plot) of annotated features with a hit
     vs annotated features without a hit */

* Junctions;

data junc_w_hit;
  set event.blast_dtct_jnc2pb_nomult_summary;
  where flag_correct_gene=1 and flag_feature_on=1;
  keep event_id;
  rename event_id=feature_id;
run;

data frags_w_hit;
   set event.pacbio_hit_correct_gene;
   keep feature_id;
run;

data feat_w_hit;
   set junc_w_hit frags_w_hit;
run;

data annot_feat;
   set event.features_w_annotations_nomulti;
   where flag_feature_on=1 and feature_type ne "fusion";
   keep feature_id feature_type flag_singleton flag_junction_annotated flag_intron_retention;
run;

proc sort data=feat_w_hit nodup;
  by feature_id;
proc sort data=annot_feat;
  by feature_id;
run;

data flag_annot_hit;
  merge annot_feat (in=in1) feat_w_hit (in=in2);
  by feature_id;
  if in2 then flag_pacbio_hit=1; else flag_pacbio_hit=0;
run;



proc freq data=flag_annot_hit;
tables flag_pacbio_hit;
run;

proc freq data=flag_annot_hit noprint;
 tables flag_pacbio_hit*flag_singleton*feature_type*flag_junction_annotated*flag_intron_retention / out=hits_by_type;
proc print data=hits_by_type;
run;

/*

                                               Cumulative
   flag_pacbio_hit    Frequency     Percent     Frequency
   --------------------------------------------------------
                 0      208042       75.46        208042
                 1       67654       24.54        275696


 flag_
pacbio_      flag_      feature_    flag_junction_    flag_intron_
  hit      singleton      type         annotated        retention     COUNT


   0           0        fragment           0                0         40264
   0           1        fragment           0                0         86532
   0           0        splicing           0                0          8208
   0           0        splicing           0                1         26688
   0           0        splicing           1                0         46350

   1           0        fragment           0                0          3796
   1           1        fragment           0                0         19388
   1           0        splicing           0                0           204
   1           0        splicing           0                1           804
   1           0        splicing           1                0         43462

*/


/* Now calc by gene */

data feat2gene;
   set event.feature2xs2gene_exp_only_nomulti;
   keep gene_id feature_id;
run;

data unannot2gene;
   set evspl.splicing_events_annot_refseq;
   keep gene_id event_id;
   rename event_id=feature_id;
run;

data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data feat2gene2;
   set unannot2gene feat2gene;
run;

proc sort data=feat2gene2 nodup;
   by gene_id feature_id;
proc sort data=exp_genes nodup;
   by gene_id;
run;

data feat2gene3;
  merge feat2gene2 (in=in1) exp_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;


proc sort data=feat2gene3 nodup;
   by feature_id gene_id;
proc sort data=flaG_annot_hit;
   by feature_id;
run;

data annot_hit_w_gene oops;
   merge flag_annot_hit (in=in1) feat2gene3 (in=in2);
   by feature_id;
   if in1 and in2 then output annot_hit_w_gene;
   else if in1 then output oops;
run;


data gene_hit_annot gene_nohit_annot gene_hit_unannot gene_nohit_unannot gene_hit_ir gene_nohit_ir;
   set annot_hit_w_gene;
   if feature_type="fragment" then do;
       if flag_pacbio_hit=1 then output gene_hit_annot;
       else output gene_nohit_annot;
       end;
   if feature_type="splicing" then do;
       if flag_junction_annotated=1 then do;
          if flag_pacbio_hit=1 then output gene_hit_annot;
          else output gene_nohit_annot;
          end;
       else do;
          if flag_intron_retention=1 then do;
             if flag_pacbio_hit=1 then output gene_hit_ir;
             else output gene_nohit_ir;
             end;
          else do;  
             if flag_pacbio_hit=1 then output gene_hit_unannot;
             else output gene_nohit_unannot;
          end;
       end;
     end;
   keep gene_id;
run;

proc sort data=gene_hit_annot nodup;
   by gene_id;
proc sort data=gene_nohit_annot nodup;
   by gene_id;
proc sort data=gene_hit_unannot nodup;
   by gene_id;
proc sort data=gene_nohit_unannot nodup;
   by gene_id;
proc sort data=gene_hit_ir nodup;
   by gene_id;
proc sort data=gene_nohit_ir nodup;
   by gene_id;
run;

data genes_w_hit;
  merge gene_hit_annot (in=in1) gene_nohit_annot (in=in2)
        gene_hit_unannot (in=in3) gene_nohit_unannot (in=in4)
        gene_hit_ir (in=in5) gene_nohit_ir (in=in6) ;
  by gene_id;
  if in1 or in3 or in5 then flag_pacbio_hit=1; else flag_pacbio_hit=0;
  if in1 or in2 then flag_annot=1; else flag_annot=0;
  if in3 or in4 then flag_unannot=1; else flag_unannot=0;
  if in5 or in6 then flag_ir=1; else flag_ir=0;
run;

proc freq data=genes_w_hit;
  tables flag_pacbio_hit;
run;

/*
                                             Cumulative    Cumulative
 flag_pacbio_hit    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0       16449       78.80         16449        78.80
               1        4426       21.20         20875       100.00

*/


proc freq data=genes_w_hit noprint;
  tables flag_pacbio_hit*flag_annot*flag_unannot*flag_ir / out=gene_count;
proc print data=gene_count;
run;


/*
  flag_
 pacbio_    flag_     flag_
   hit      annot    unannot    flag_ir    COUNT

    0         1         0          0       10989
    0         1         0          1        3859
    0         1         1          0         322
    0         1         1          1        1279
    1         1         0          0         445
    1         1         0          1        1458
    1         1         1          0         135
    1         1         1          1        2388

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
  by flag_pacbio_hit feature_type flag_junction_annotated flag_intron_retention;
proc means data=annot_hit_w_apn noprint;
  by flag_pacbio_hit feature_type  flag_junction_annotated flag_intron_retention;
  var mean_apn_npc;
  output out=apn_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc print data=apn_distrib ; run;

/*

       flag_
      pacbio_  feature_  flag_junction_  flag_intron_
 Obs    hit      type       annotated      retention   _TYPE_  _FREQ_          mean            sd

  1      0     fragment         0              0          0    126796  33.186075398  263.14217635
  2      0     splicing         0              0          0      8208  0.7144877391  2.4776635254
  3      0     splicing         0              1          0     26688  6.4648316413  243.63292119
  4      0     splicing         1              0          0     46350  15.104570254  122.50466866

  5      1     fragment         0              0          0     23184  104.35429151  353.93787073
  6      1     splicing         0              0          0       204  3.0076602598  10.320013166
  7      1     splicing         0              1          0       804  7.1697380918  27.654912965
  8      1     splicing         1              0          0     43462  55.995473699  277.60633293



 Obs           min              q1          median              q3             max

  1   0.0001643115    0.3414634146     1.762245102    10.200922266    25968.547727
  2    0.024691358    0.2407407407    0.3395061728    0.5679012346    134.17283951
  3   0.0182926829    0.3231707317    0.3414634146    0.9756097561    24728.530488
  4    0.024691358    0.3456790123    1.0185185185    2.8456790124     4569.345679
  5   0.0138888889    11.453003381    31.142351915    85.490197568     13274.41954
  6   0.1419753086    0.3796296296    1.0216049383    2.0648148148    136.79012346
  7   0.2012195122    0.9207317073    2.3414634146    5.9176829268    503.76829268
  8   0.1728395062    4.1419753086    12.231481481    37.037037037    19137.549383



Okay, looks like the APN is generally higher for features with PB hits.

The unannotated events (unannotated junctions, IR events) aren't too different though, maybe marginally higher
*/

/* Make permenant and export for plots */

data event.pacbio_hits_feat_mean_apn;
   set annot_hit_w_apn;
   drop _TYPE_ _FREQ_;
run;


proc export data=event.pacbio_hits_feat_mean_apn
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_hits_features_apn.csv"
     dbms=csv replace;
run;


