/* Fusions summary by gene */

/* set libraries */
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Get list of genes with num_exons_total num_exons_exp num_DE */

data genes_de_summary;
  set sugrue.results_gene_summary;
  keep gene_id exon_count total_exon_count de_exon_cnt updown_cat;
  rename updown_cat=exon_updown_cat;
run;

/* get fusions and exons differentially expressed and cat to gene level */

* get differentially expressed fusions first;

data fusions_for_genes;
   set sugrue.results_by_fusions_final;
   if flag_multigene=0 and flag_p05=1;
   keep fusion_id gene_id;
run;

proc sort data=fusions_for_genes;
  by gene_id;
run;

* Cat fusions and exons by gene;
*	Gene, Flag_1_DE_Fusion, Flag_1_DE_SE, Flag_1_DE_IR;
*DE_fusion >= DE_splicing >= DE_IR;
* Nominal P<0.05 as few FDR-corrected exons significant;

* get counts first;

proc freq noprint data=fusions_for_genes;
   tables gene_id / out=fus4gene_cnt;
run;

proc sort data=fus4gene_cnt;
  by descending count;
run;
*max=47 fusions per gene;

proc sort data=fusions_for_genes;
   by gene_cat fusion_id;
run;

data fusions_for_genes2; 
  array fusions[47] $ 14;
  array exons[47] $ 3234;
  retain fusions1-fusions47;
  retain exons1-exons47;

  set fusions_for_genes;
  by gene_cat;
  
  if first.gene_cat then do;
     call missing(of fusions1-fusions47);
     call missing(of exons1-exons47);
     records = 0;
  end;

  records + 1;
  fusions[records]=fusion_id;
  exons[records]=exon_cat;
  if last.gene_cat then output;
run;

  *clean up the output file;

data fusions_for_genes3;
  set fusions_for_genes2;
  length fusions_cat $ 705;
  length exons_cat $ 2241;
  fusions_cat= cats(OF fusions1-fusions47);
  exons_cat= cats(OF exons1-exons47);
  drop records fusions1-fusions47 fusion_id exon_cat exons1-exons47 ;
  rename gene_cat=gene_id;
  run;

* sort and merge with gene summary;

proc sort data=genes_de_summary;
   by gene_id;
run;

proc sort data=fusions_for_genes3;
   by gene_id;
run;

data genes_de_summary_w_fus oops_no_exp;
   merge genes_de_summary (in=in1) fusions_for_genes3 (in=in2);
   by gene_id;
   if in1 and in2 then output genes_de_summary_w_fus;
   else if in1 then do;
       fusions_cat='.';
       exons_cat='.';
       output genes_de_summary_w_fus;
       end;
   else output oops_no_exp;
run;

/* Make permenant for now */

data sugrue.genes_de_summary_w_fus;
   set genes_de_summary_w_fus;
run;

