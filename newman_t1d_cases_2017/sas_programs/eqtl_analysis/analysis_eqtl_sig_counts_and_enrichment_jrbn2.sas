
/* Get counts and test for T1D enrichment */


libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* get counts for FDR 5% */

data eqtl_summary;
  set eqtl.eqtl_results_summary_clean;
  keep snp_id feature_id feature_type flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
  if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
    if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
      if flag_cd19_fdr05=. then flag_cd19_fdr05=0;
  run;

proc sort data=eqtl_summary nodup;
    by snp_id feature_id;
    run;
        
proc freq data=eqtl_summary noprint;
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_all;
   run;
   
proc freq data=eqtl_summary noprint;
   where feature_type='exon';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_exon;
   run;
      
   proc freq data=eqtl_summary noprint;
   where feature_type='AS';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_as;
   run;
   
   proc freq data=eqtl_summary noprint;
   where feature_type='IR';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_ir;
   run;
   

/* SUMMARY
	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	ALL
ALL	1906	1936	2063	333	326	219	964
EXON	840	751	740	174	151	111	425
JUNC	967	1060	1183	146	160	105	507
IR	99	125	140	13	15	3	32

NEW COUNTS:
	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	ALL
ALL	1929	1977	2080	359	337	219	970
EXON	861	788	756	200	160	110	431
JUNC	968	1064	1184	146	162	106	507
IR	100	125	140	13	15	3	32

*/

/* Get counts on gene - FDR 5% */

data eqtl_summary_by_gene_fus eqtl_summary_by_gene_splice; 
  set eqtl.eqtl_results_summary_clean;

  if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
  if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
  if flag_cd19_fdr05=. then flag_cd19_fdr05=0;

  keep gene_id feature_id snp_id feature_type flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;

  if feature_type='exon' then output eqtl_summary_by_gene_fus;
  else output eqtl_summary_by_gene_splice;
    run;

data fusion2gene;
   set fus.unique_info_fusions_si;
   rename fusion_id=feature_id;
run;

data eqtl_summary_by_gene_fus2;
   set eqtl_summary_by_gene_fus;
   drop gene_id;
run;

proc sort data=fusion2gene;
   by feature_id;
proc sort data=eqtl_summary_by_gene_fus2;
   by feature_id;
run;


data eqtl_summary_by_gene_fus3;
   merge eqtl_summary_by_gene_fus2 (in=in1) fusion2gene (in=in2);
   by feature_id;
   if in1;
run;

data eqtl_summary_by_gene;
   set eqtl_summary_by_gene_fus3 eqtl_summary_by_gene_splice;
run;

proc sort data=eqtl_summary_by_gene out=eqtl_summary_by_gene2 nodup;
   by gene_id feature_id snp_id ;
   run;

proc means data=eqtl_summary_by_gene2 noprint;
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
   output out=eqtls_summed_by_gene_all sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='exon';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
   output out=eqtls_summed_by_gene_exon sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='AS';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
   output out=eqtls_summed_by_gene_as sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='IR';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
   output out=eqtls_summed_by_gene_ir sum=;
run;


data eqtls_summed_by_gene_all2;
   set eqtls_summed_by_gene_all;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_exon2;
   set eqtls_summed_by_gene_exon;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_as2;
   set eqtls_summed_by_gene_as;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_ir2;
   set eqtls_summed_by_gene_ir;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;


/* Merge in T1D gene flags */

data immunoflags;
   set con.immunogene_flags;
   keep gene_id flag_diabetes_gene;
   run;

proc sort data=eqtls_summed_by_gene_all2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_exon2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_as2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_ir2;
   by gene_id;
proc sort data=immunoflags;
   by gene_id;
   run;



data eqtls_by_gene_all;
   merge eqtls_summed_by_gene_all2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;
   
data eqtls_by_gene_exon;
   merge eqtls_summed_by_gene_exon2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;

data eqtls_by_gene_as;
   merge eqtls_summed_by_gene_as2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;

data eqtls_by_gene_ir;
   merge eqtls_summed_by_gene_ir2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;
   
   
/* Get counts for each subset */

proc freq data=eqtls_by_gene_all noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_all;
   run;

proc freq data=eqtls_by_gene_exon noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_exon;
   run;
   
proc freq data=eqtls_by_gene_as noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_as;
   run;

proc freq data=eqtls_by_gene_ir noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_ir;
   run;

/* SUMMARY

OLD:
	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	ALL
ALL	90	89	101	56	52	46	172
EXON	71	82	82	43	33	22	91
JUNC	54	55	71	39	43	35	97
IR	26	38	30	7	9	4	14
					
NEW:
	CD4	CD8	CD19	CD4/8	CD4/19	CD8/19	ALL
ALL	95	94	111	59	54	48	176
EXON	75	87	92	46	36	23	95
JUNC	54	55	70	38	43	36	98
IR	27	38	30	7	9	4	14

*/


