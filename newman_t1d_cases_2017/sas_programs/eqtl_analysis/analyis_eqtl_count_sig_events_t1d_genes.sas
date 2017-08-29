libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';


/* Get genes with a sig eQTL */


data exon4 exon8 exon19 junc4 junc8 junc19 ir4 ir8 ir19 all4 all8 all19;
   set eqtl.results_summary_table_w_means;
   /* CD4 */
   if CD4_FDR_P ne . then do;
         output all4;
         if CD4_FDR_P lt 0.05 and feature_type='exon' then output exon4;
         if CD4_FDR_P lt 0.05 and feature_type='Junc' then output junc4;
         if CD4_FDR_P lt 0.05 and feature_type='IR' then output ir4;
         end;

   if CD8_FDR_P ne . then do;
         output all8;
         if CD8_FDR_P lt 0.05 and feature_type='exon' then output exon8;
         if CD8_FDR_P lt 0.05 and feature_type='Junc' then output junc8;
         if CD8_FDR_P lt 0.05 and feature_type='IR' then output ir8;
         end;

   if CD19_FDR_P ne . then do;
         output all19;
         if CD19_FDR_P lt 0.05 and feature_type='exon' then output exon19;
         if CD19_FDR_P lt 0.05 and feature_type='Junc' then output junc19;
         if CD19_FDR_P lt 0.05 and feature_type='IR' then output ir19;
         end;

   keep gene_id;
run;

data eqtl_gene;
   set eqtl.results_summary_table_w_means;
   if CD19_FDR_P ne . or CD4_FDR_P ne . or  CD8_FDR_P ne . ;
   keep gene_id;
run;


data ibase;
  set con.immunobase_gene_flags;
  if flag_immunobase_diabetes_gene=1;
  keep gene_id;
run;

%macro count_eqtl(dataset,header);

   proc freq data=&dataset. noprint;
      tables gene_id / out=&dataset._2;
   run;

   data &dataset._3;
     set &dataset._2;
     keep gene_id count;
     rename count=&header.;
   run;

   proc sort data=&dataset._3;
   by gene_id;
   run;
%mend;


%count_eqtl(all4,all_cd4);
%count_eqtl(exon4,exon_cd4);
%count_eqtl(junc4,junc_cd4);
%count_eqtl(ir4,ir_cd4);

%count_eqtl(all8,all_cd8);
%count_eqtl(exon8,exon_cd8);
%count_eqtl(junc8,junc_cd8);
%count_eqtl(ir8,ir_cd8);

%count_eqtl(all19,all_cd19);
%count_eqtl(exon19,exon_cd19);
%count_eqtl(junc19,junc_cd19);
%count_eqtl(ir19,ir_cd19);

proc sort data=ibase nodup;
  by gene_id;
proc sort data=eqtl_gene nodup;
  by gene_id;
run;

data eqtl_counts_t1d;
   merge ibase (in=in1) eqtl_gene (in=in2) all4_3 all8_3 all19_3 
         exon4_3 exon8_3 exon19_3
         junc4_3 junc8_3 junc19_3
         ir4_3 ir8_3 ir19_3;
         by gene_id;
         if in1 and in2;
run;

/* replace missing with 0 */

data eqtl_counts_t1d2;
   set eqtl_counts_t1d;
   array change _numeric_;
      do over change;
      if change=. then change=0;
            end;
   run ;


/* Make permenant */

data eqtl.t1d_eqtl_table_summary_v2;
   set eqtl_counts_t1d2;
run;

/* Export */

proc export data=eqtl.t1d_eqtl_table_summary_v2 outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/t1d_eqtls_summary_table_v3.csv'
   dbms=csv replace;
run; quit;

/* Count */

data count;
   set eqtl.t1d_eqtl_table_summary_v2;
   if exon_cd4 gt 0 or exon_cd8 gt 0 or exon_cd19 gt 0 then flag_exon=1; else flag_exon=0;
   if junc_cd4 gt 0 or junc_cd8 gt 0 or junc_cd19 gt 0 then flag_junc=1; else flag_junc=0;
   if ir_cd4 gt 0 or ir_cd8 gt 0 or ir_cd19 gt 0 then flag_ir=1; else flag_ir=0;
run;

proc freq data=count;
   tables flag_exon flag_junc flag_ir flag_exon*flag_junc flag_exon*flag_ir flag_junc*flag_ir flag_exon*flag_junc*flag_ir;
run;


/* check overlap between new eQTL genes and revised Ibase list */

data eqtl_gene;
 set eqtl.results_w_onengut_credible_means;
 if num_credible_snps_80perc gt 0 and flag_eqtl_sig=1;
 keep gene_id;
run;

*unstack;


data eqtl_gene2;
   length gene_id2 $35.; 
   set eqtl_gene;
   do i=1 by 1 while(scan(gene_id,i,'|') ^=' '); 
gene_id2=scan(gene_id,i,'|'); 
keep gene_id2;
output; 
end; 
rename gene_id2=gene_id;
run;


data ibase;
  set con.immunobase_gene_flags;
  if flag_immunobase_diabetes_gene=1;
  keep gene_id;
run;

proc sort data=eqtl_gene2 nodup;
  by gene_id;
proc sort data=ibase nodup;
  by gene_id;
run;

data eqtl ibase eqtl_ibase;
  merge ibase (in=in1) eqtl_gene2 (in=in2);
  by gene_id;
  if in1 and in2 then output eqtl_ibase;
  else if in1 then output ibase;
  else output eqtl;
run;

