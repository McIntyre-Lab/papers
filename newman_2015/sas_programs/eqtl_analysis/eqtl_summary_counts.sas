/* Counts for eQTLs */
/* Need:
       total SNPs selected
       number of genes with expression
       expressed features
               exons
               junctions
               IR
	Combinations:
		Total
		Exons	
		Junctions
		IR
 */


/* Total selected SNPs */

data snps;
  set eqtl.eqtl_results_summary_table;
  keep snp_id;
run;

proc sort data=snps nodup;
  by snp_id;
run;


/* Total genes */

data genes;
  set eqtl.eqtl_results_summary_table;
  keep gene_id;
run;

data genes2;
length gene_id2 $36.; 
   set genes;
   do i=1 by 1 while(scan(gene_id,i,'|') ^=' ');
      gene_id2=scan(gene_id,i,'|');
      output;
   end;
   keep gene_id2;
   rename gene_id2=gene_id;
run;

proc sort data=genes2 nodup;
  by gene_id;
run;

data immunoflags;
   set con.immunogene_flags;
   if flag_immuno_gene=1;
   keep gene_id;
run;

proc sort data=immunoflags;
   by gene_id;
run;

data immunogenes;
  merge immunoflags (in=in1) genes2 (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Total features */

data features;
  set eqtl.eqtl_results_summary_table;
  keep feature_id;
run;

proc sort data=features nodup;
  by feature_id;
run;


/* Total exons */

data exons;
  set eqtl.eqtl_results_summary_table;
  if feature_type='exon';
  keep feature_id;
run;

proc sort data=exons nodup;
  by feature_id;
run;

/* Total junctions */

data juncs;
  set eqtl.eqtl_results_summary_table;
  if feature_type='Junc';
  keep feature_id;
run;

proc sort data=juncs nodup;
  by feature_id;
run;


/* Total IR */
data ir;
  set eqtl.eqtl_results_summary_table;
  if feature_type='IR';
  keep feature_id;
run;

proc sort data=ir nodup;
  by feature_id;
run;



/* All-3 exons */

data exons;
  set eqtl.eqtl_results_summary_table;
  if feature_type='exon' and CD4_FDR_P ne . and CD8_FDR_P ne . and CD19_FDR_P ne .;
  keep feature_id;
run;

proc sort data=exons nodup;
  by feature_id;
run;

/* All-3 junctions */

data juncs;
  set eqtl.eqtl_results_summary_table;
  if feature_type='Junc' and CD4_FDR_P ne . and CD8_FDR_P ne . and CD19_FDR_P ne .;
  keep feature_id;
run;

proc sort data=juncs nodup;
  by feature_id;
run;


/* All-3 IR */
data ir;
  set eqtl.eqtl_results_summary_table;
  if feature_type='IR' and CD4_FDR_P ne . and CD8_FDR_P ne . and CD19_FDR_P ne .;
  keep feature_id;
run;

proc sort data=ir nodup;
  by feature_id;
run;

/* Total combos */

data combos_all;
   set eqtl.eqtl_results_summary_table;
   keep feature_id snp_id;
run;

proc sort data=combos_all nodup;
  by feature_id snp_id;
run;


/* Exon combos */

data combos_exon;
   set eqtl.eqtl_results_summary_table;
   if feature_type='exon';
   keep feature_id snp_id;
run;

proc sort data=combos_exon nodup;
  by feature_id snp_id;
run;



/* Junction combos */

data combos_junc;
   set eqtl.eqtl_results_summary_table;
   if feature_type='Junc';
   keep feature_id snp_id;
run;

proc sort data=combos_junc nodup;
  by feature_id snp_id;
run;


/* IR combos */

data combos_ir;
   set eqtl.eqtl_results_summary_table;
   if feature_type='IR';
   keep feature_id snp_id;
run;

proc sort data=combos_ir nodup;
  by feature_id snp_id;
run;


/* Num sig, num genes */

data num_sig_cd4;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P lt 0.05;
   cell_type='CD4';
   keep feature_id snp_id gene_id cell_type;
run;

data num_sig_cd8;
   set eqtl.eqtl_results_summary_table;
   if CD8_FDR_P lt 0.05;
   cell_type='CD8';
   keep feature_id snp_id gene_id cell_type;
run;

data num_sig_cd19;
   set eqtl.eqtl_results_summary_table;
   if CD19_FDR_P lt 0.05;
   cell_type='CD19';
   keep feature_id snp_id gene_id cell_type;
run;

data num_sig_tests;
   set num_sig_cd4 num_sig_cd8 num_sig_cd19;
run;

data num_sig_count;
  set eqtl.eqtl_results_summary_table;
  if flag_eqtl_sig=1;
run;


/* eQTL counts for events only detected in all tissues */

data eqtl_clean_data;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P=. then delete;
   if CD8_FDR_P=. then delete;
   if CD19_FDR_P=. then delete;
   if CD4_FDR_P lt 0.05 then flag_cd4_sig=1; else flag_cd4_sig=0;
   if CD8_FDR_P lt 0.05 then flag_cd8_sig=1; else flag_cd8_sig=0;
   if CD19_FDR_P lt 0.05 then flag_cd19_sig=1; else flag_cd19_sig=0;
run;

proc freq data=eqtl_clean_data noprint;
   tables flag_cd4_sig*flag_cd8_sig*flag_cd19_sig / out=all_dtct_all;
run;

proc freq data=eqtl_clean_data noprint;
   where feature_type='exon';
   tables flag_cd4_sig*flag_cd8_sig*flag_cd19_sig / out=all_dtct_exon;
run;

proc freq data=eqtl_clean_data noprint;
   where feature_type='Junc';
   tables flag_cd4_sig*flag_cd8_sig*flag_cd19_sig / out=all_dtct_junc;
run;

proc freq data=eqtl_clean_data noprint;
   where feature_type='IR';
   tables flag_cd4_sig*flag_cd8_sig*flag_cd19_sig / out=all_dtct_ir;
run;


/* Single-tissue specifc features with associations */
 
data feat_single_cell;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P ne . and CD8_FDR_P eq . and CD19_FDR_P eq . and flag_eqtl_sig=1 then output;
   if CD4_FDR_P eq . and CD8_FDR_P ne . and CD19_FDR_P eq . and flag_eqtl_sig=1 then output;
   if CD4_FDR_P eq . and CD8_FDR_P eq . and CD19_FDR_P ne . and flag_eqtl_sig=1 then output;
run;


/* CD4 */

data feat_cd4_only;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P ne . and CD8_FDR_P eq . and CD19_FDR_P eq . and flag_eqtl_sig=1 then output;
run;


/* CD8 */
data feat_cd8_only;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P eq . and CD8_FDR_P ne . and CD19_FDR_P eq . and flag_eqtl_sig=1 then output;
run;


/* CD19 */
data feat_cd19_only;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P eq . and CD8_FDR_P eq . and CD19_FDR_P ne . and flag_eqtl_sig=1 then output;
run;



/* Two-tissue specifc features with associations in both */

data feat_two_tissue_both;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P ne . and CD8_FDR_P ne . and CD19_FDR_P eq . then do;
           if CD4_FDR_P lt 0.05 and CD8_FDR_P lt 0.05 then output; end;
   if CD4_FDR_P ne . and CD8_FDR_P eq . and CD19_FDR_P ne . then do;
           if CD4_FDR_P lt 0.05 and CD19_FDR_P lt 0.05 then output; end;
   if CD4_FDR_P eq . and CD8_FDR_P ne . and CD19_FDR_P ne . then do;
           if CD8_FDR_P lt 0.05 and CD19_FDR_P lt 0.05 then output; end;
run;



/* Two-tissue specifc features with associations in one */

data feat_two_tissue_one;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P ne . and CD8_FDR_P ne . and CD19_FDR_P eq . then do;
           if CD4_FDR_P lt 0.05 and CD8_FDR_P ge 0.05 then output;
           if CD4_FDR_P ge 0.05 and CD8_FDR_P lt 0.05 then output;
           end;
   if CD4_FDR_P ne . and CD8_FDR_P eq . and CD19_FDR_P ne . then do;
           if CD4_FDR_P lt 0.05 and CD19_FDR_P ge 0.05 then output;
           if CD4_FDR_P ge 0.05 and CD19_FDR_P lt 0.05 then output;
           end;
   if CD4_FDR_P eq . and CD8_FDR_P ne . and CD19_FDR_P ne . then do;
           if CD8_FDR_P lt 0.05 and CD19_FDR_P ge 0.05 then output;
           if CD8_FDR_P ge 0.05 and CD19_FDR_P lt 0.05 then output;
           end;
run;


/* Remaining # */
 
data feat_all_cell;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P eq . then delete;
   if CD8_FDR_P eq . then delete;
   if CD19_FDR_P eq . then delete;
   if flag_eqtl_sig=1 then output;
run;

/* Total -- to check */

data feat_all_cell;
   set eqtl.eqtl_results_summary_table;
   if flag_eqtl_sig=1 then output;
run;



/* Genes with all */

