
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Determine the set of genes that are detected. Then for the set of genes detected in all three cell types,
   count the number "on". The number that are DD are now evidence of alternative splicing!

   I would like to count the following:

   cd4*cd8*cd19 for fusions/events from genes on in all cell types

   cd4*cd8 for fusions/events from genes only on in CD4 and CD8
   cd4*cd19 for fusions/events from genes only on in CD4 and CD19
   cd8*cd19 for fusions/events from genes only on in CD8 and CD19

   cd4 fusions/events from genes only on in CD4
   cd8 fusions/events from genes only on in CD8
   cd19 fusions/events from genes only on in CD19

Then we can look at fusions/events from genes detected in at least 2 cell type

*/

data fusion_exp;
  set con.fusions_on_gt_apn0;
  keep fusion_id flag_cd19_on flag_cd4_on flag_cd8_on;
run;

*Since not dealing with multigene splicing events, need to summarize fusions info to individual genes, and not multigenes;
data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=fus2gene;
   by fusion_id;
proc sort data=fusion_exp;
   by fusion_id;
run;

data fus2gene_w_exp oops;
   merge fus2gene (in=in1) fusion_exp (in=in2);
   by fusion_id;
   if in1 and in2 then output fus2gene_w_exp;
   else output oops;
run;

*Summarize expression flag on gene;

proc sort data=fus2gene_w_exp;
   by gene_id;
run;

proc means data=fus2gene_w_exp noprint;
   by gene_id;
   var flag_cd19_on flag_cd4_on flag_cd8_on;
   output out=gene_exp
          sum(flag_cd4_on)=num_exons_exp_cd4
          sum(flag_cd8_on)=num_exons_exp_cd8
          sum(flag_cd19_on)=num_exons_exp_cd19;
run;

/* Flag if gene is on */

data flag_gene_on;
  set gene_exp;
  if num_exons_exp_cd4 > 0 then flag_cd4_gene_on=1; else flag_cd4_gene_on=0;
  if num_exons_exp_cd8 > 0 then flag_cd8_gene_on=1; else flag_cd8_gene_on=0;
  if num_exons_exp_cd19 > 0 then flag_cd19_gene_on=1; else flag_cd19_gene_on=0;
  drop _TYPE_ _FREQ_;
run;

/* Count gene-level flags */

proc freq data=flag_gene_on noprint;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on / out=gene_count;
proc print data=gene_count;
run;

/*
                            flag_
 flag_cd4_    flag_cd8_     cd19_
  gene_on      gene_on     gene_on    COUNT

     0            0           0       25321
     0            0           1        1804
     0            1           0         626
     0            1           1         467
     1            0           0         648
     1            0           1         497
     1            1           0        1994
     1            1           1       41027

*/

/* Make permenant */

data con.flag_gene_detection_by_cell;
   set flag_gene_on;
run;

