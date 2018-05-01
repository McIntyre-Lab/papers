ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";
libname con "!PATCON/sas_data";
libname hg19 "!PATCON/useful_human_data/aceview_hg19/fusions/sas_data";


/* For each gene (nomulti!!), count the number of exons that are detected in both cell types at APN>0
   and the number of exons that are detected in both cell types are APN>5

   use fusions for this!! */


data fus_on_apn0;
   set con.fusions_on_gt_apn0;
   /* Count fusions that are on in each combination of cell types */
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=1 then flag_fusion_all_on=1; else flag_fusion_all_on=0;

   if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_cd4_on=1; else flag_fusion_cd4_on=0;
   if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=0 then flag_fusion_cd8_on=1; else flag_fusion_cd8_on=0;
   if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=1 then flag_fusion_cd19_on=1; else flag_fusion_cd19_on=0;

   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=0 then flag_fusion_cd4cd8_on=1; else flag_fusion_cd4cd8_on=0;
   if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=1 then flag_fusion_cd4cd19_on=1; else flag_fusion_cd4cd19_on=0;
   if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=1 then flag_fusion_cd8cd19_on=1; else flag_fusion_cd8cd19_on=0;
   
   keep fusion_id flag_fusion_all_on flag_fusion_cd4_on flag_fusion_cd8_on flag_fusion_cd19_on
                  flag_fusion_cd4cd8_on flag_fusion_cd4cd19_on flag_fusion_cd8cd19_on ;
run;

proc freq data=fus_on_apn0 noprint;
   tables flag_fusion_all_on*flag_fusion_cd4_on*flag_fusion_cd8_on*flag_fusion_cd19_on*
                  flag_fusion_cd4cd8_on*flag_fusion_cd4cd19_on*flag_fusion_cd8cd19_on / out=gene_check;
run;

proc print data=gene_check;
run; *okay, this looks correct;


data fus2gene;
  set hg19.hg19_aceview_fusions_si_info;
  keep fusion_id gene_id;
run;

data gene_w_multi;
  set hg19.hg19_aceview_fusions_si_info;
  where flag_multigene=1;
  keep gene_id;
run;


proc sort data=fus2gene nodup;
  by gene_id fusion_id;
proc sort data=gene_w_multi nodup;
  by gene_id;
run;

data fus2gene_nomulti;
  merge fus2gene (in=in1) gene_w_multi (in=in2);
  by gene_id;
  if in2 then delete;
  else if in1 then output;
run;


proc sort data=fus_on_apn0;
   by fusion_id;
proc sort data=fus2gene_nomulti nodup;
  by fusion_id gene_id;
run;

data fus2gene_on;
  merge fus2gene_nomulti (in=in1) fus_on_apn0 (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2gene_on;
  by gene_id;
proc means data=fus2gene_on noprint;
  by gene_id;
  var flag_fusion_all_on flag_fusion_cd4_on flag_fusion_cd8_on flag_fusion_cd19_on
                  flag_fusion_cd4cd8_on flag_fusion_cd4cd19_on flag_fusion_cd8cd19_on ;
  output out=gene_on sum=;
run;

data flag_gene_on;
   set gene_on;
   /* Now flag genes that are on in each cell types. The end should be 23236 genes exp in all */

   if sum(flag_fusion_all_on,flag_fusion_cd4_on,flag_fusion_cd4cd8_on,flag_fusion_cd4cd19_on) > 0 
      then flag_gene_cd4_on=1; else flag_gene_cd4_on=0;

   if sum(flag_fusion_all_on,flag_fusion_cd8_on,flag_fusion_cd4cd8_on,flag_fusion_cd8cd19_on) > 0 
      then flag_gene_cd8_on=1; else flag_gene_cd8_on=0;

   if sum(flag_fusion_all_on,flag_fusion_cd19_on,flag_fusion_cd4cd19_on,flag_fusion_cd8cd19_on) > 0 
      then flag_gene_cd19_on=1; else flag_gene_cd19_on=0;

run;

proc freq data=flag_gene_on noprint;
   tables flag_gene_cd4_on*flag_gene_cd8_on*flag_gene_cd19_on / out=gene_count;
run;

proc print data=gene_count;
run;



/*
 flag_     flag_     flag_
 gene_     gene_     gene_
cd4_on    cd8_on    cd19_on    COUNT

   0         0         0       21721
   0         0         1        1442
   0         1         0         511
   0         1         1         360
   1         0         0         502
   1         0         1         379
   1         1         0        1563
   1         1         1       23236

*/

/* Make permenant */

data event.t1d_flag_gene_on_by_cell;
  set flag_gene_on;
run;


