/* Make a summary table of just the T1D candidates */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';


/* Want:

T1D genes only

Summary of total signifcant
CD4
CD8
CD19

Stratify:
 exon, junction, ir
 - split on cell type

*/

data t1d_eqtls;
   set eqtl.eqtl_results_summary_table;
   if cd4_fdr_p=. then flag_cd4_sig=.; else if cd4_fdr_p lt 0.05 then flag_cd4_sig=1; else flag_cd4_sig=0;
   if cd8_fdr_p=. then flag_cd8_sig=.; else if cd8_fdr_p lt 0.05 then flag_cd8_sig=1; else flag_cd8_sig=0;
   if cd19_fdr_p=. then flag_cd19_sig=.; else if cd19_fdr_p lt 0.05 then flag_cd19_sig=1; else flag_cd19_sig=0;
   if flag_diabetes_gene=1 then output;
run;

/* Stack gene ids */

data t1d_eqtls_stacked(rename=new=gene_id2);
length new $35.; 
set t1d_eqtls; 
do i=1 by 1 while(scan(gene_id,i,'|') ^=' '); 
new=scan(gene_id,i,'|'); 
output; 
end; 
run;

data t1d_eqtls_stacked2;
   set t1d_eqtls_stacked;
   drop i gene_id flag_diabetes_gene flag_autoimmune_gene flag_multigene;
   rename gene_id2=gene_id;
run;

data t1d_gene_list;
  set con.immunogene_flags;
  if flag_diabetes_gene=1;
  keep gene_id;
run;

proc sort data=t1d_eqtls_stacked2;
   by gene_id;
proc sort data=t1d_gene_list;
   by gene_id;
run;

data t1d_gene_eqtls not_t1d;
   merge t1d_eqtls_stacked2 (in=in1) t1d_gene_list (in=in2);
   by gene_id;
   if in1 and in2 then output t1d_gene_eqtls;
   else if in1 then output not_t1d;
run;

proc freq data=t1d_gene_eqtls noprint;
   tables gene_id*flag_eqtl_sig / out=t1d_sig_eqtls;
   tables gene_id*flag_cd4_sig / out=t1d_sig_eqtls_cd4;
   tables gene_id*flag_cd8_sig / out=t1d_sig_eqtls_cd8;
   tables gene_id*flag_cd19_sig / out=t1d_sig_eqtls_cd19;
run;


proc freq data=t1d_gene_eqtls noprint;
   where feature_type='exon';
   tables gene_id*flag_eqtl_sig / out=t1d_sig_eqtls_exon;
   tables gene_id*flag_cd4_sig / out=t1d_sig_eqtls_cd4_exon;
   tables gene_id*flag_cd8_sig / out=t1d_sig_eqtls_cd8_exon;
   tables gene_id*flag_cd19_sig / out=t1d_sig_eqtls_cd19_exon;
run;

proc freq data=t1d_gene_eqtls noprint;
   where feature_type='Junc';
   tables gene_id*flag_eqtl_sig / out=t1d_sig_eqtls_junc;
   tables gene_id*flag_cd4_sig / out=t1d_sig_eqtls_cd4_junc;
   tables gene_id*flag_cd8_sig / out=t1d_sig_eqtls_cd8_junc;
   tables gene_id*flag_cd19_sig / out=t1d_sig_eqtls_cd19_junc;
run;

proc freq data=t1d_gene_eqtls noprint;
   where feature_type='IR';
   tables gene_id*flag_eqtl_sig / out=t1d_sig_eqtls_ir;
   tables gene_id*flag_cd4_sig / out=t1d_sig_eqtls_cd4_ir;
   tables gene_id*flag_cd8_sig / out=t1d_sig_eqtls_cd8_ir;
   tables gene_id*flag_cd19_sig / out=t1d_sig_eqtls_cd19_ir;
run;

%macro format_tables(dataset,cell_type,feature);

data &dataset._2;
  set &dataset.;
  if flag_&cell_type._sig=1;
  keep gene_id count;
  rename count=count_&cell_type._sig_&feature.;
run;

proc sort data=&dataset._2;
   by gene_id;
run;

%mend;

%format_tables(t1d_sig_eqtls,eqtl,all);
%format_tables(t1d_sig_eqtls_cd4,cd4,all);
%format_tables(t1d_sig_eqtls_cd8,cd8,all);
%format_tables(t1d_sig_eqtls_cd19,cd19,all);

%format_tables(t1d_sig_eqtls_exon,eqtl,exon);
%format_tables(t1d_sig_eqtls_cd4_exon,cd4,exon);
%format_tables(t1d_sig_eqtls_cd8_exon,cd8,exon);
%format_tables(t1d_sig_eqtls_cd19_exon,cd19,exon);

%format_tables(t1d_sig_eqtls_junc,eqtl,junc);
%format_tables(t1d_sig_eqtls_cd4_junc,cd4,junc);
%format_tables(t1d_sig_eqtls_cd8_junc,cd8,junc);
%format_tables(t1d_sig_eqtls_cd19_junc,cd19,junc);

%format_tables(t1d_sig_eqtls_ir,eqtl,ir);
%format_tables(t1d_sig_eqtls_cd4_ir,cd4,ir);
%format_tables(t1d_sig_eqtls_cd8_ir,cd8,ir);
%format_tables(t1d_sig_eqtls_cd19_ir,cd19,ir);


data t1d_eqtl_summary_all;
   merge t1d_sig_eqtls_2 t1d_sig_eqtls_cd4_2 t1d_sig_eqtls_cd8_2 t1d_sig_eqtls_cd19_2;
   by gene_id;
run;

data t1d_eqtl_summary_exon;
   merge t1d_sig_eqtls_exon_2 t1d_sig_eqtls_cd4_exon_2 t1d_sig_eqtls_cd8_exon_2 t1d_sig_eqtls_cd19_exon_2;
   by gene_id;
run;


data t1d_eqtl_summary_junc;
   merge  t1d_sig_eqtls_junc_2 t1d_sig_eqtls_cd4_junc_2 t1d_sig_eqtls_cd8_junc_2 t1d_sig_eqtls_cd19_junc_2;
   by gene_id;
run;


data t1d_eqtl_summary_ir;
   merge  t1d_sig_eqtls_ir_2 t1d_sig_eqtls_cd4_ir_2 t1d_sig_eqtls_cd8_ir_2 t1d_sig_eqtls_cd19_ir_2;
   by gene_id;
run;

data t1d_eqtl_summary;
   merge t1d_eqtl_summary_all t1d_eqtl_summary_exon t1d_eqtl_summary_junc t1d_eqtl_summary_ir;
   by gene_id;
run;
 
data t1d_eqtl_summary2;
   set t1d_eqtl_summary;
   array change _numeric_;
      do over change;
      if change=. then change=0;
            end;
   run ;


/* Make permenant */

data eqtl.t1d_eqtl_table_summary;
   set t1d_eqtl_summary2;
run;


/* Export */

proc export data=eqtl.t1d_eqtl_table_summary outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/t1d_eqtls_summary_table_v2.csv'
   dbms=csv replace;
run; quit;


