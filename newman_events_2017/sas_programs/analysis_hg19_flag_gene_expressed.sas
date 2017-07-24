ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* For the case-only data, flag if gene is expressed and if it has multigene regions
   Here, I am saying a gene is "expressed" if it has at least one fusion in at least one cell type
   detected at APN>0
*/

data gene_w_multi;
   set hg19.hg19_aceview_fusions_si_info;
   where flag_multigene=1;
   keep gene_id;
run;

data gene_exp;
   set con.flag_gene_detection_by_cell;
   if flag_cd4_gene_on=1 or flag_cd8_gene_on=1 or flag_cd19_gene_on=1;
   keep gene_id flag_cd4_gene_on flag_cd8_gene_on flag_cd19_gene_on;
run;

proc sort data=gene_w_multi nodup;
  by gene_id;
proc sort data=gene_exp nodup;
  by gene_id;
run;

data gene_flags;
   merge gene_exp (in=in1) gene_w_multi (in=in2);
   by gene_id;
   if in2 then flag_multigene=1; else flag_multigene=0;
   if in1 then output;
run;

proc freq data=gene_flags noprint;
   tables flag_multigene*flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on / out=gene_count;
proc print data=gene_count;
run;

/*
                                           flag_
     flag_      flag_cd4_    flag_cd8_     cd19_
   multigene     gene_on      gene_on     gene_on    COUNT    PERCENT

These have cell-specific expression.
Keep these so we can check that transcript estimates make sense
(ie, these are "control genes", in that we expect these to only have
estimates for one cell type):
       0            0            0           1        1442     3.0640
       0            0            1           0         511     1.0858
       0            1            0           0         502     1.0667

Keep these:
       0            0            1           1         360     0.7649
       0            1            0           1         379     0.8053
       0            1            1           0        1563     3.3211
       0            1            1           1       23236    49.3721

But not these:
       1            0            0           1         362     0.7692
       1            0            1           0         115     0.2444
       1            0            1           1         107     0.2274
       1            1            0           0         146     0.3102
       1            1            0           1         118     0.2507
       1            1            1           0         431     0.9158
       1            1            1           1       17791    37.8025
*/

data event.hg19_flag_gene_expressed;
   set gene_flags;
run;


