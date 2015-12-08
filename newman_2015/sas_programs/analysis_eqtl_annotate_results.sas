/* Make eQTL summaries */

libname eqtl '/mnt/data/eqtls/sas_data';
libname splice '/mnt/data/splice';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';



data splicing_annot;
   set splice.splicing_events_annotations;
   length feature_type $4.;
   if flag_intron_retention=1 then feature_type='IR';
   else feature_type='AS';
   keep event_id feature_type gene_id;
   rename event_id=feature_id;
run;

proc sort data=eqtl.eqtl_results_w_fdr;
   by feature_id;
proc sort data=splicing_annot;
   by feature_id;
run;


data eqtl_splicing_events eqtl_exon_events;
   merge eqtl.eqtl_results_w_fdr (in=in1) splicing_annot (in=in2);
   by feature_id;
   if in1 and in2 then output eqtl_splicing_events;
   else if in1 then do;
	feature_type="exon";
        output eqtl_exon_events;
        end;
run;


data fusion2gene;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id;
   rename fusion_id=feature_id;
run;

proc sort data=fusion2gene;
   by feature_id;
proc sort data=eqtl_exon_events;
   by feature_id;
run;

data eqtl_exon_events2;
   set eqtl_exon_events;
   drop gene_id;
run;

data eqtl_exon_events3;
  merge eqtl_exon_events2 (in=in1) fusion2gene (in=in2);
  by feature_id;
  if in1;
run;


data eqtl_results_w_genes;
   set eqtl_splicing_events eqtl_exon_events2;
run;



data cd4_test cd8_test cd19_test cd4_sig cd8_sig cd19_sig;
   set eqtl_results_w_genes;
   if cell_type='CD4' then output cd4_test;
   else if cell_type='CD8' then output cd8_test;
   else if cell_type='CD19' then output cd19_test;

   if cell_type='CD4' and flag_eqtl_fdr05=1 then output cd4_sig;
   else if cell_type='CD8' and flag_eqtl_fdr05=1 then output cd8_sig;
   else if cell_type='CD19' and flag_eqtl_fdr05=1 then output cd19_sig;
run;

data eqtl_test_index;
   set eqtl_results_w_genes;
   keep feature_id snp_id gene_id feature_type;
run;

proc sort data=eqtl_test_index nodup;
   by feature_id snp_id cell_type;
run;


data cd4_test2;
   set cd4_test;
   flag_cd4_tested=1;
   keep feature_id snp_id ProbF fdr_p flag_cd4_tested;
   rename ProbF=CD4_pvalue fdr_p=cd4_fdr;
run;

data cd8_test2;
   set cd8_test;
   flag_cd8_tested=1;
   keep feature_id snp_id ProbF fdr_p flag_cd8_tested;
   rename ProbF=CD8_pvalue fdr_p=cd8_fdr;
run;


data cd19_test2;
   set cd19_test;
   flag_cd19_tested=1;
   keep feature_id snp_id ProbF fdr_p flag_cd19_tested;
   rename ProbF=CD19_pvalue fdr_p=cd19_fdr;
run;

data cd4_sig2;
   set cd4_sig;
   flag_cd4_fdr05=1;
   keep feature_id snp_id flag_cd4_fdr05;
   run;

data cd8_sig2;
   set cd8_sig;
   flag_cd8_fdr05=1;
   keep feature_id snp_id flag_cd8_fdr05;
   run;


data cd19_sig2;
   set cd19_sig;
   flag_cd19_fdr05=1;
   keep feature_id snp_id flag_cd19_fdr05;
   run;

proc sort data=eqtl_test_index;
   by feature_id snp_id;
proc sort data=cd4_test2;
   by feature_id snp_id;
proc sort data=cd8_test2;
   by feature_id snp_id;
proc sort data=cd19_test2;
   by feature_id snp_id;
proc sort data=cd4_sig2;
   by feature_id snp_id;
proc sort data=cd8_sig2;
   by feature_id snp_id;
proc sort data=cd19_sig2;
   by feature_id snp_id;
run;


data eqtl_results_summary;
  merge eqtl_test_index cd4_test2 (in=in1) cd8_test2 (in=in2) cd19_test2 (in=in3) cd4_sig2 (in=in4) cd8_sig2 (in=in5) cd19_sig2 (in=in6);
   by feature_id snp_id;
run;



data eqtl_results_summary2;
   set eqtl_results_summary;
    if flag_cd4_tested=. then flag_cd4_tested=0;
    if flag_cd8_tested=. then flag_cd8_tested=0;
    if flag_cd19_tested=. then flag_cd19_tested=0;
    if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
    if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
    if flag_cd19_fdr05=. then flag_cd19_fdr05=0;
run;

data eqtl.eqtl_results_summary;
   set eqtl_results_summary2;
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=0 and flag_cd19_fdr05=0 then flag_cell_specific=1;
   else if flag_cd4_fdr05=0 and flag_cd8_fdr05=1 and flag_cd19_fdr05=0 then flag_cell_specific=1;
   else if flag_cd4_fdr05=0 and flag_cd8_fdr05=0 and flag_cd19_fdr05=1 then flag_cell_specific=1;
   else flag_cell_specific=0;
run;

