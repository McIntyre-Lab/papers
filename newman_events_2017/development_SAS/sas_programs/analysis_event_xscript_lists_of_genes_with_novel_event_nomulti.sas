ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* From my list of transcripts with at least one unique feature detected OR have no assigned unique features,
   I want to count the number of these transcripts that correspond to:
     genes with at least one novel junction
     genes with at least one novel IR

   Count the number of transcripts, and the number of genes */


data exp_genes;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

data genes_novel;
   set event.unannotated_events_by_gene;
run;

proc sort data=exp_genes;
  by gene_id;
proc sort data=genes_novel;
  by gene_id;
run;

data exp_genes_w_novel;
  merge exp_genes (in=in1) genes_novel (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc freq data=exp_genes_w_novel;
  tables flag_gene_has_unannotated_junc*flag_gene_has_ir;
run;

*457 genes with unannotated junctions;
*5317 genes with IR;
*3667 genes with IR and unannotated junctions;

data xs2gene;
  set event.feature2xs2gene_exp_only_nomulti;
  keep transcript_id gene_id;
run;

proc sort data=exp_genes_w_novel;
   by gene_id;
proc sort data=xs2gene nodup;
   by gene_id transcript_id;
run;

data exp_novel_w_xs;
   merge xs2gene (in=in1) exp_genes_w_novel (in=in2);
   by gene_id;
   if in1 and in2;
run;

data prob_xs;
   set event.feature_dtct_cnt_by_xscript_exp;
   where bin_xscript_perc_uniq_dtct ne "0%";
   keep transcript_id;
run;

proc sort data=prob_xs;
   by transcript_id;
proc sort data=exp_novel_w_xs;
   by transcript_id;
run;

data exp_novel_w_prob_xs;
   merge exp_novel_w_xs (in=in1) prob_xs (in=in2);
   by transcript_id;
   if in1 and in2;
run;

proc freq data=exp_novel_w_prob_xs;
    tables flag_gene_has_unannotated_junc flag_gene_has_ir;
run;


/*
   flag_gene_has_                             Cumulative    Cumulative
 unannotated_junc    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       17020       50.02         17020        50.02
                1       17008       49.98         34028       100.00


                                              Cumulative    Cumulative
 flag_gene_has_ir    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1668        4.90          1668         4.90
                1       32360       95.10         34028       100.00



*/

proc freq data=exp_genes_w_novel;
  tables flag_gene_has_unannotated_junc flag_gene_has_ir;
run;


/*
    flag_gene_has_                             Cumulative    Cumulative
  unannotated_junc    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        5317       56.32          5317        56.32
                 1        4124       43.68          9441       100.00


                                               Cumulative    Cumulative
  flag_gene_has_ir    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         457        4.84           457         4.84
                 1        8984       95.16          9441       100.00


*/


