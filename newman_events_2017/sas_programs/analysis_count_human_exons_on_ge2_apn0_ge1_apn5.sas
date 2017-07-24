ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname hg19 "!PATCON/useful_human_data/aceview_hg19/fusions/sas_data";
libname con "!PATCON/sas_data";

/* For each gene (nomulti!!), count the number of exons that are detected in both cell types at APN>0
   and the number of exons that are detected in both cell types are APN>5

   use fusions for this!! */

/* Flag genes with multigene fusions */

data genes_W_multi;
  set hg19.hg19_aceview_fusions_si_info;
  where flag_multigene=1;
  keep gene_id;
run;

data fus_on_apn0;
   set con.fusions_on_gt_apn0;
   keep fusion_id flag_cd19_on flag_cd8_on flag_cd4_on flag_fusion_all_on0;
run;


data fus_on_apn5;
   set con.fusion_on_ge_apn5_v2;
   if flag_cd4_on=. or flag_cd8_on=. or flag_cd19_on=. then delete;
   keep fusion_id flag_cd19_on flag_cd8_on flag_cd4_on flag_fusion_all_on5;
   rename flag_cd19_on=flag_cd19_on_apn5
          flag_cd8_on=flag_cd8_on_apn5
          flag_cd4_on=flag_cd4_on_apn5;
run;

data fus2gene;
  set hg19.hg19_aceview_fusions_si_info;
  where flag_multigene=0;
  keep gene_id fusion_id;
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
   var flag_cd19_on flag_cd8_on flag_cd4_on flag_fusion_all_on0
       flag_cd19_on_apn5 flag_cd8_on_apn5 flag_cd4_on_apn5 flag_fusion_all_on5;
   output out=fus_on_by_gene sum=;
run;

data flag_genes;
   set fus_on_by_gene;
   if flag_cd19_on ge 2 then flag_cd19_ge2_apn0=1; else flag_cd19_ge2_apn0=0;
   if flag_cd8_on ge 2  then flag_cd8_ge2_apn0=1; else flag_cd8_ge2_apn0=0;
   if flag_cd4_on ge 2  then flag_cd4_ge2_apn0=1; else flag_cd4_ge2_apn0=0;
   if flag_fusion_all_on0 ge 2  then flag_all_ge2_apn0=1; else flag_all_ge2_apn0=0;
   if flag_cd19_on_apn5 ge 1 then flag_cd19_ge1_apn5=1; else flag_cd19_ge1_apn5=0;
   if flag_cd8_on_apn5 ge 1 then flag_cd8_ge1_apn5=1; else flag_cd8_ge1_apn5=0;
   if flag_cd4_on_apn5 ge 1 then flag_cd4_ge1_apn5=1; else flag_cd4_ge1_apn5=0;
   if flag_fusion_all_on5 ge 1 then flag_all_ge1_apn5=1; else flag_all_ge1_apn5=0;
run;

proc sort data=flag_genes;
   by gene_id;
proc sort data=genes_w_multi nodup;
   by gene_id;
run;

data flag_genes2;
  merge flag_genes (in=in1) genes_w_multi (in=in2);
   by gene_id;
   if in2 then flag_multigene=1; else flag_multigene=0;
   if in1;
run;

proc freq data=flag_genes2 noprint;
  tables flag_cd19_ge2_apn0*flag_cd4_ge2_apn0*flag_cd8_ge2_apn0*flag_all_ge2_apn0*
         flag_cd19_ge1_apn5*flag_cd4_ge1_apn5*flag_cd8_ge1_apn5*flag_all_ge1_apn5*flag_multigene / out=gene_count;
run;
proc print data=gene_count;
run;

/* Save genes and fusions to keep for further analysis */


data gene2keep_all;
   set flag_genes2;
   where flag_all_ge1_apn5=1 and flaG_all_ge2_apn0=1;
   keep gene_id;
run;

data gene2keep_cd4cd8;
   set flag_genes2;
   where flag_cd4_ge1_apn5=1 and flaG_cd4_ge2_apn0=1
     and flag_cd8_ge1_apn5=1 and flaG_cd8_ge2_apn0=1;
   keep gene_id;
run;

data gene2keep_cd4cd19;
   set flag_genes2;
   where flag_cd4_ge1_apn5=1 and flaG_cd4_ge2_apn0=1
     and flag_cd19_ge1_apn5=1 and flaG_cd19_ge2_apn0=1;
   keep gene_id;
run;

data gene2keep_cd8cd19;
   set flag_genes2;
   where flag_cd19_ge1_apn5=1 and flaG_cd19_ge2_apn0=1
     and flag_cd8_ge1_apn5=1 and flaG_cd8_ge2_apn0=1;
   keep gene_id;
run;


data fus2keep_all;
   set fus_on2gene;
   where flag_fusion_all_on0=1 ;
   keep gene_id fusion_id;
run;

data fus2keep_cd4cd8;
   set fus_on2gene;
   where flag_cd4_on=1 and flag_cd8_on=1 ;
   keep gene_id fusion_id;
run;

data fus2keep_cd4cd19;
   set fus_on2gene;
   where flag_cd4_on=1 and flag_cd19_on=1 ;
   keep gene_id fusion_id;
run;

data fus2keep_cd8cd19;
   set fus_on2gene;
   where flag_cd19_on=1 and flag_cd8_on=1 ;
   keep gene_id fusion_id;
run;

proc sort data=gene2keep_all;
  by gene_id;
proc sort data=fus2keep_all;
  by gene_id;
proc sort data=gene2keep_cd4cd8;
  by gene_id;
proc sort data=fus2keep_cd4cd8;
  by gene_id;
proc sort data=gene2keep_cd4cd19;
  by gene_id;
proc sort data=fus2keep_cd4cd19;
  by gene_id;
proc sort data=gene2keep_cd8cd19;
  by gene_id;
proc sort data=fus2keep_cd8cd19;
  by gene_id;
run;

data gene_fus2keep_all;
  merge gene2keep_all (in=in1) fus2keep_all (in=in2);
  by gene_id;
  if in1 and in2;
run;

data gene_fus2keep_cd4cd8;
  merge gene2keep_cd4cd8 (in=in1) fus2keep_cd4cd8 (in=in2);
  by gene_id;
  if in1 and in2;
run;

data gene_fus2keep_cd4cd19;
  merge gene2keep_cd4cd19 (in=in1) fus2keep_cd4cd19 (in=in2);
  by gene_id;
  if in1 and in2;
run;


data gene_fus2keep_cd8cd19;
  merge gene2keep_cd8cd19 (in=in1) fus2keep_cd8cd19 (in=in2);
  by gene_id;
  if in1 and in2;
run;

data event.hg19_genes_for_ds_test_all;
   set gene_fus2keep_all;
run;

data event.hg19_genes_for_ds_test_cd4cd8;
   set gene_fus2keep_cd4cd8;
run;

data event.hg19_genes_for_ds_test_cd4cd19;
   set gene_fus2keep_cd4cd19;
run;

data event.hg19_genes_for_ds_test_cd8cd19;
   set gene_fus2keep_cd8cd19;
run;


