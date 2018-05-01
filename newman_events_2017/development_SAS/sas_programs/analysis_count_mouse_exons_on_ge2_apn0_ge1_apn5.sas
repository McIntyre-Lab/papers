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

data fus_on_apn5;
   set event.flag_fusion_on_apn5;
   if flag_fusion_nsc_on=. then delete;
   if flag_fusion_old_on=. then delete;
   if flag_fusion_nsc_on=1 and flag_fusion_old_on=1 then flag_fusion_both_on_apn5=1;
   else flaG_fusion_both_on_apn5=0;
   rename flag_fusion_nsc_on=flag_fusion_nsc_on_apn5
          flag_fusion_old_on=flag_fusion_old_on_apn5;
run;

data fus2gene;
  set event.feature2xs2gene_nomulti;
  keep feature_id gene_id;
  rename feature_id=fusion_id;
run;

proc sort data=fus_on_apn0;
   by fusion_id;
proc sort data=fus_on_apn5;
   by fusion_id;
proc sort data=fus2gene nodup;
   by fusion_id gene_id;
run;

data fus_on2gene;
  merge fus2gene (in=in1) fus_on_apn0 (in=in2) fus_on_apn5;
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus_on2gene;
   by gene_id;
proc means data=fus_on2gene noprint;
   by gene_id;
   var flag_fusion_nsc_on flag_fusion_old_on flag_fusion_both_on
       flag_fusion_nsc_on_apn5 flag_fusion_old_on_apn5 flag_fusion_both_on_apn5;
   output out=fus_on_by_gene
sum(flag_fusion_nsc_on)=nsc_on_apn0
sum(flag_fusion_old_on)=old_on_apn0
sum(flag_fusion_both_on)=both_on_apn0
sum(flag_fusion_nsc_on_apn5)=nsc_on_apn5
sum(flag_fusion_old_on_apn5)=old_on_apn5
sum(flag_fusion_both_on_apn5)=both_on_apn5;
run;

data flag_genes;
   set fus_on_by_gene;
   if nsc_on_apn0 ge 2 then flag_nsc_ge2_apn0=1; else flag_nsc_ge2_apn0=0; 
   if old_on_apn0 ge 2 then flag_old_ge2_apn0=1; else flag_old_ge2_apn0=0; 
   if both_on_apn0 ge 2 then flag_both_ge2_apn0=1; else flag_both_ge2_apn0=0; 

   if nsc_on_apn5 ge 1 then flag_nsc_ge1_apn5=1; else flag_nsc_ge1_apn5=0; 
   if old_on_apn5 ge 1 then flag_old_ge1_apn5=1; else flag_old_ge1_apn5=0; 
   if both_on_apn5 ge 1 then flag_both_ge1_apn5=1; else flag_both_ge1_apn5=0; 
run;

proc freq data=flag_genes noprint;
  tables flaG_nsc_ge2_apn0*flaG_old_ge2_apn0*flaG_both_ge2_apn0*
         flaG_nsc_ge1_apn5*flaG_old_ge1_apn5*flaG_both_ge1_apn5 / out=gene_count;
run;
proc print data=gene_count;
run;


/*
  flag_       flag_       flag_       flag_       flag_       flag_
nsc_ge2_    old_ge2_      both_     nsc_ge1_    old_ge1_      both_
  apn0        apn0      ge2_apn0      apn5        apn5      ge1_apn5    COUNT

    0           0           0           0           0           0       11745
    0           0           0           0           1           0          13
    0           0           0           1           0           0          66
    0           0           0           1           1           1         374
    0           1           0           0           0           0        1353
    0           1           0           0           1           0          39
    0           1           0           1           1           1           6
    1           0           0           0           0           0        1764
    1           0           0           0           1           0           1
    1           0           0           1           0           0           2
    1           0           0           1           1           1          11
    1           1           0           0           0           0         422
    1           1           0           1           1           1           1

    1           1           1           0           0           0        4971
    1           1           1           0           1           0         801
    1           1           1           1           0           0         625

    1           1           1           1           1           1        6272 <- these ones!

*/

/* Save these genes and fusions for further analysis */

data gene2keep;
   set flag_genes;
   where flag_both_ge1_apn5=1 and flaG_both_ge2_apn0=1;
   keep gene_id;
run;

data fus2keep;
   set fus_on2gene;
   where flag_fusion_both_on=1 ;
   keep gene_id fusion_id;
run;

proc sort data=gene2keep;
  by gene_id;
proc sort data=fus2keep;
  by gene_id;
run;

data gene_fus2keep;
  merge gene2keep (in=in1) fus2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

data event.mm10_genes_for_ds_test;
   set gene_fus2keep;
run;


