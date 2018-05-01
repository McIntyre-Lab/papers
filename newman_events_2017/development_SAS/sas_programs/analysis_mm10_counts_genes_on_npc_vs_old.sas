ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";

/* For each gene (nomulti!!), count the number of exons that are detected in both cell types at APN>0
   and the number of exons that are detected in both cell types are APN>5

   use fusions for this!! */


data fus_on_apn0;
   set event.flag_fusion_on;
   if flag_fusion_nsc_on=1 and flag_fusion_old_on=1 then flag_fusion_both_on=1;
   else flaG_fusion_both_on=0;
run;

data fus2gene;
  set event.feature2xs2gene_nomulti;
  keep feature_id gene_id;
  rename feature_id=fusion_id;
run;

proc sort data=fus_on_apn0;
   by fusion_id;
proc sort data=fus2gene nodup;
  by fusion_id gene_id;
run;

data fus2gene_on;
  merge fus2gene (in=in1) fus_on_apn0 (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2gene_on;
  by gene_id;
proc means data=fus2gene_on noprint;
  by gene_id;
  var flag_fusion_nsc_on flag_fusion_old_on flag_fusion_both_on ;
  output out=gene_on sum=;
run;

data flag_gene_on;
   set gene_on;
   if flag_fusion_nsc_on>0 then flag_gene_nsc_on=1; else flag_gene_nsc_on=0;
   if flag_fusion_old_on>0 then flag_gene_old_on=1; else flag_gene_old_on=0;
   if flag_gene_nsc_on=1 and flag_gene_old_on=1 then flag_gene_both_on=1; else flag_gene_both_on=0;
run;

proc freq data=flag_gene_on;
   tables flag_gene_nsc_on*flag_gene_old_on flag_gene_both_on;
run;


/*

  flag_fusion_nsc_on
            flag_fusion_old_on

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   6195 |   1396 |   7591
           |  21.76 |   4.90 |  26.67
           |  81.61 |  18.39 |
           |  68.33 |   7.20 |
  ---------+--------+--------+
         1 |   2871 |  18004 |  20875
           |  10.09 |  63.25 |  73.33
           |  13.75 |  86.25 |
           |  31.67 |  92.80 |
  ---------+--------+--------+
  Total        9066    19400    28466
              31.85    68.15   100.00


      flag_gene_                             Cumulative    Cumulative
         both_on    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       10462       36.75         10462        36.75
               1       18004       63.25         28466       100.00


*/

/* Make permenant */

data event.mm10_flag_gene_on_by_cell;
  set flag_gene_on;
run;


