libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Update diabetes gene flag then re-export */

data eqtl_genes;
    set eqtl.results_summary_table_w_means;
    keep gene_id;
run;


proc sort data=eqtl_genes nodup;
   by gene_id;
run;


* stack IDs;
data eqtl_gene2;
   length gene_id2 $35.; 
   set eqtl_genes;
   do i=1 by 1 while(scan(gene_id,i,'|') ^=' '); 
gene_id2=scan(gene_id,i,'|'); 
drop i;
output; 
end; 
run;


/* Get ibase flags */

data ibase;
   set con.immunobase_gene_flags;
   if flag_immunobase_diabetes_gene=1;
   keep gene_id;
   rename gene_id=gene_id2;
run;

proc sort data=ibase;
   by gene_id2;
proc sort data=eqtl_gene2;
   by gene_id2;
run;

data eqtl2ibase;
   merge eqtl_gene2 (in=in1) ibase (in=in2);
   by gene_id2;
   if in2 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;

data diabetes_genes;
   set eqtl2ibase;
   if flag_immunobase_diabetes_gene=1;
   keep gene_id;
run;

proc sort data=diabetes_genes nodup;
  by gene_id;
proc sort data=eqtl.results_summary_table_w_means;
  by gene_id;
proc sort data=eqtl.results_w_onengut_credible_means;
  by gene_id;
run;

data results_summary_table_w_means;
   merge diabetes_genes (in=in1) eqtl.results_summary_table_w_means (in=in2);
   by gene_id;
   if in1 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in2 then output;
run;

data results_w_onengut_credible_means;
   merge diabetes_genes (in=in1) eqtl.results_w_onengut_credible_means (in=in2);
   by gene_id;
   if in1 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in2 then output;
run;

data eqtl.results_summary_table_w_means_v2;
   retain gene_id flag_multigene flag_autoimmune_gene flag_immunobase_diabetes_gene;
   set results_summary_table_w_means;
   drop flag_diabetes_gene;
   rename flag_immunobase_diabetes_gene=flag_diabetes_gene;
run;

data eqtl.results_w_onengut_cred_means_v2;
   retain gene_id snp_id onengut_index_snp_cat flag_immunobase_diabetes_gene;
   set results_w_onengut_credible_means;
   drop flag_diabetes_gene;
   rename flag_immunobase_diabetes_gene=flag_diabetes_gene;
run;

proc export data=eqtl.results_summary_table_w_means_v2
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/AJHG_submission_Aug2016/SuppTable1.csv' dbms=csv replace;
run;

proc export data=eqtl.results_w_onengut_cred_means_v2
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/AJHG_submission_Aug2016/SuppTable3.csv' dbms=csv replace;
run;



