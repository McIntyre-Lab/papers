ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For BLAST hits of unannotated features, I want to count the number of total detected unannotated
   features with hits (and by unannotated junction/IR)
   and the number of genes with unannotated hits (and by unannotated junction/IR) */

data features_w_hit;
   set event.pacbio_hit_correct_gene;
   keep feature_id;
run;

data unannot_feat;
   set event.features_w_annotations_nomulti;
   where flag_feature_on=1 and feature_type="splicing";
   keep feature_id flag_intron_retention flag_junction_annotated;
run;

proc sort data=features_w_hit nodup;
  by feature_id;
proc sort data=unannot_feat;
  by feature_id;
run;

data unannot_w_hit;
  merge  unannot_feat (in=in1)  features_w_hit (in=in2);
  by feature_id;
  if in2 then flag_pacbio_hit=1; else flag_pacbio_hit=0;
  if in1 then output;
run;

proc freq data=unannot_w_hit;
  tables flag_junction_annotated*flag_pacbio_hit;
run;


/*
   flag_junction_annotated
             flag_pacbio_hit

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  35904 |      0 |  35904
            |  28.56 |   0.00 |  28.56
            | 100.00 |   0.00 |
            |  43.13 |   0.00 |
   ---------+--------+--------+
          1 |  47337 |  42475 |  89812
            |  37.65 |  33.79 |  71.44
            |  52.71 |  47.29 |
            |  56.87 | 100.00 |
   ---------+--------+--------+
   Total       83241    42475   125716
               66.21    33.79   100.00

No hits for unannotated events!
*/

/* Now I want to plot the distribution of unannotated junctions/IR (do separately) and compare then to
   annotated junctions with/without PB hits */

data event_counts;
    set event.mm10_refseq_splicing_counts;
    where sample_id ? "NSC";
    keep sample_id event_id apn;
    rename event_id=feature_id;
run;

proc sort data=event_counts;
   by feature_id sample_id;
proc means data=event_counts noprint;
   by feature_id;
   var apn;
   output out=mean_apn_by_feat mean=mean_apn_npc;
run;

data unannot_w_hit2;
  set unannot_w_hit;
  if flag_junction_annotated=1 then delete;
run;


proc sort data=mean_apn_by_feat;
  by feature_id;
proc sort data=unannot_w_hit2;
   by feature_id;
run;

data unannot_hit_w_apn;
  merge unannot_w_hit2 (in=in1) mean_apn_by_feat (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* before I export, I want to check the distribution of APN here */

proc sort data=unannot_hit_w_apn;
  by flag_pacbio_hit flag_intron_retention;
proc means data=unannot_hit_w_apn noprint;
  by flag_pacbio_hit flag_intron_retention;
  var mean_apn_npc;
  output out=apn_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc print data=apn_distrib ; run;

/*
        flag_
       pacbio_    flag_intron_
Obs      hit        retention     _TYPE_    _FREQ_            mean              sd

 1        0             0            0        8412    0.7700996261    2.9470025863
 2        0             1            0       27492    6.4854465397    240.09040531



Obs              q1          median              q3             max

 1     0.2407407407    0.3395061728    0.5987654321    136.79012346
 2     0.3292682927    0.4024390244    1.0121951219    24728.530488


Okay, looks like the APN is generally higher for features with PB hits
*/

/* Export for plots */


proc export data=unannot_hit_w_apn
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_hits_unannotated_features_apn.csv"
     dbms=csv replace;
run;


