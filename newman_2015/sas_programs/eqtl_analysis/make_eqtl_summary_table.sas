/* Make a summary table of eQTL results */

libname con '/home/jrbnewman/concannon/sas_data';
libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';


/* We want columns:
gene_id flag_mutligene (for fusions) flag_autoimmune_gene flag_diabetes_gene feature_type feature_id snp_id flag_eqtl_sig CD4_FDR_P CD8_FDR_P CD19_FDR_P

Then sort by putting T1D gene eqtls first and export. Put T1D IR events first */

/* Trim results summary */

data results_summary_trim;
   set eqtl.eqtl_results_summary_clean;
   if flag_cd4_fdr05=1 or flag_cd8_fdr05=1 or flag_cd19_fdr05=1 then flag_eqtl_sig=1;
   else flag_eqtl_sig=0;
   keep gene_id feature_id snp_id feature_type flag_eqtl_sig cd4_fdr cd8_fdr cd19_fdr;
run;

/* Add gene_id back into exon eQTLs */

data exon_eqtls splicing_eqtls;
   set results_summary_trim;
   if feature_type='exon' then output exon_eqtls;
   else output splicing_eqtls;
run;

/* Add flag_immunogenes for individual fusions */


data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

data immunogenes;
   set con.immunogene_flags;
   if flag_immuno_gene=1;
   keep gene_id flag_immuno_gene;
run;

data diabetes_genes;
   set con.immunogene_flags;
   if flag_diabetes_gene=1;
   keep gene_id flag_diabetes_gene;
run;


proc sort data=fus2gene nodup;
   by gene_id;
proc sort data=immunogenes;
   by gene_id;
proc sort data=diabetes_genes;
   by gene_id;
run;


data immuno_fusions;
   merge fus2gene (in=in1) immunogenes (in=in2);
   by gene_id;
   if in1 and in2;
run;

data t1d_fusions;
   merge fus2gene (in=in1) diabetes_genes (in=in2);
   by gene_id;
   if in1 and in2;
run;

data immuno_fusions2;
   set immuno_fusions;
   keep fusion_id flag_immuno_gene;
run;

data t1d_fusions2;
   set t1d_fusions;
   keep fusion_id flag_diabetes_gene;
run;

proc sort data=immuno_fusions2 nodup;
   by fusion_id;
proc sort data=t1d_fusions2 nodup;
   by fusion_id;
run;


data fus_w_immunoflags oops;
   merge immuno_fusions2 (in=in1) t1d_fusions2 (in=in2);
   by fusion_id;
   if in1 and in2 then output fus_w_immunoflags;
   else if in1 then do;
       flag_diabetes_gene=1;
       output fus_w_immunoflags;
       end;
   else output oops;
run;

/* Get fusions and flag multigene */

data fus2multigene;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id flag_multigene;
run;

proc sort data=fus2multigene;
   by fusion_id;
proc sort data=fus_w_immunoflags;
   by fusion_id;
run;


data fus2gene_w_flags oops;
   merge fus2multigene (in=in1) fus_w_immunoflags (in=in2);
   by fusion_id;
   if in1 and in2 then output fus2gene_w_flags;
   else if in1 then do;
       flag_immuno_gene=0;
       flag_diabetes_gene=0;
       output fus2gene_w_flags;
       end;
   else output oops;
run;


/* Merge flags in with exons eqtls */

data fusions_tested;
   set exon_eqtls;
   keep feature_id;
run;

data fus2gene_w_flags2;
   length fusion_id $2550.;
   format fusion_id $2475.;
   set fus2gene_w_flags;
   rename fusion_id=feature_id;
run;

proc sort data=fusions_tested nodup;
   by feature_id;
proc sort data=fus2gene_w_flags2;
   by feature_id;
run;

data tested_exons_w_flags no_eqtl oops;
   merge fusions_tested (in=in1) fus2gene_w_flags2 (in=in2);
   by feature_id;
   if in1 and in2 then output tested_exons_w_flags;
   else if in1 then output oops;
   else output no_eqtl;
run;

data exon_eqtls2;
   set exon_eqtls;
   drop gene_id;
run;

proc sort data=tested_exons_w_flags;
   by feature_id;
proc sort data=exon_eqtls2;
   by feature_id;
run;

data exon_eqtls_w_flags oops1 oops2;
   merge exon_eqtls2 (in=in1) tested_exons_w_flags (in=in2);
   by feature_id;
   if in1 and in2 then output exon_eqtls_w_flags;
   else if in1 then output oops1;
   else output oops2;
run;


/* Merge immunogene flags with splicing eqtls */

data immunogene_flags;
   set con.immunogene_flags;
   keep gene_id flag_immuno_gene flag_diabetes_gene;
run;

proc sort data=immunogene_flags;
   by gene_id;
proc sort data=splicing_eqtls;
   by gene_id;
run;

data splicing_eqtls_w_flags no_eqtls oops;
   merge splicing_eqtls (in=in1) immunogene_flags (in=in2);
   by gene_id;
   if in1 and in2 then output splicing_eqtls_w_flags;
   else if in1 then output oops;
   else output no_eqtls;
run;

/* Stack eQTL results */

data eqtls_w_flags;
   set exon_eqtls_w_flags splicing_eqtls_w_flags (in=in2);
   if in2 then flag_multigene=0;
run;


/* Remove duplicates */

proc sort data=eqtls_w_flags out=eqtls_w_flags_nodup nodup;
   by feature_id snp_id gene_id;
run;


/* Make permenant */

data eqtl.eqtl_results_summary_table;
   retain gene_id flag_multigene flag_immuno_gene flag_diabetes_gene feature_type feature_id snp_id flag_eqtl_sig CD4_FDR CD8_FDR CD19_FDR;
   set eqtls_w_flags_nodup;
   if feature_type='AS' then feature_type='Junc';
   rename flag_immuno_gene=flag_autoimmune_gene cd4_fdr=CD4_FDR_P cd8_fdr=CD8_FDR_P cd19_fdr=CD19_FDR_P;
run;

proc sort data=eqtl.eqtl_results_summary_table;
   by descending flag_eqtl_sig descending flag_diabetes_gene descending flag_autoimmune_gene;
run;

/* Export */

proc export data=eqtl.eqtl_results_summary_table outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_summary_table.csv'
    dbms=csv replace;
run; quit;

proc export data=eqtl.eqtl_results_summary_table outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/eqtl_summary_table.csv'
    dbms=csv replace;
run; quit;

proc sort data=eqtl.eqtl_results_summary_table;
   by descending flag_eqtl_sig descending flag_diabetes_gene descending flag_autoimmune_gene feature_type;
run;


proc export data=eqtl.eqtl_results_summary_table outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_summary_table_ir_first.csv'
    dbms=csv replace;
run; quit;

proc export data=eqtl.eqtl_results_summary_table outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/eqtl_summary_table_ir_first.csv'
    dbms=csv replace;
run; quit;
