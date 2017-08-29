libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splice '/mnt/store/splice';
libname splicing '/mnt/data/splicing';
ods listing; ods html close;
libname splictmp '/mnt/store/diabetes_sandbox/splicing';


data lsmeans;
  set splicing.anova_lsmeans;
  keep event_id cell_type estimate;
run;

proc sort data=lsmeans;
  by event_id cell_type;
proc transpose data=lsmeans out=lsmeans_sbys;
  by event_id;
  var estimate;
  id cell_type;
run;

data splicing.anova_lsmeans_diffs;
  set lsmeans_sbys;
  if CD4 > CD8 then magnitude_cd4cd8=-CD4/CD8; else magnitude_cd4cd8=CD8/CD4;
  if CD4 > CD19 then magnitude_cd4cd19=-CD4/CD19; else magnitude_cd4cd19=CD19/CD4;
  if CD8 > CD19 then magnitude_cd8cd19=-CD8/CD19; else magnitude_cd8cd19=CD19/CD8;

  keep event_id magnitude_cd4cd8 magnitude_cd4cd19 magnitude_cd8cd19 CD4 CD8 CD19;
  rename CD4=lsmeans_CD4 CD8=lsmeans_CD8 CD19=lsmeans_CD19;
run;

/* Full exon results table:

chr fusion_start fusion_stop fusion_id gene_id exon_id flag_multigene
flag_immuno_gene flag_ibase_diabetes
flag_cd19_on flag_cd4_on flag_cd8_onflag_all_on
mean_apn_cd4 (w/ SD) mean_apn_cd8 (w/ SD) mean_apn_cd19 (w/ SD)
contrast_cd4cd8_p contrast_cd4cd8_fdr
contrast_cd4cd19_p contrast_cd4cd19_fdr
contrast_cd8cd19_p contrast_cd8cd19_fdr */

/* Get splicing annotations */

data splice_info;
  set splice.splicing_events_annotations;
  if feature1_stop-feature1_start > 37 then feature1_start=feature1_stop - 37;
  if feature2_stop-feature2_start > 37 then feature2_stop=feature2_start + 37;
  if num_transcripts > 0 then flag_junction_annotated=1;
  keep event_id gene_id chr strand num_transcripts transcript_id feature1_id 
  feature1_type feature2_id feature2_type feature1_start feature1_stop
  feature2_start feature2_stop flag_junction_annotated flag_intron_retention
  flag_exonskip flag_alt_donor flag_alt_acceptor;
run;


/* Convert gene-level flags to fusion-level flags */

data ai_genes t1d_genes;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai_genes;
   if flag_immunobase_diabetes_gene=1 then output t1d_genes;
   keep gene_id;
run;

proc sort data=splice_info nodup;
  by gene_id;
proc sort data=ai_genes nodup;
  by gene_id;
proc sort data=t1d_genes nodup;
  by gene_id;
run;

data splice_info_w_flags;
  merge splice_info (in=in1) ai_genes (in=in2) t1d_genes (in=in3);
  by gene_id;
  if in2 then flag_autoimmune_gene=1; else flag_autoimmune_gene=0;
  if in3 then flag_diabetes_gene=1; else flag_diabetes_gene=0;
  if in1 then output;
run;


/* Add in "on" flags */

data on_flags;
  set splicing.splicing_results_clean;
  keep event_id flag_cd19_on flag_cd4_on flag_cd8_on;
run;

proc sort data=on_flags;
  by event_id;
proc sort data=splice_info_w_flags;
  by event_id;
run;

data splice_info_w_flags2 oops;
  merge splice_info_w_flags (in=in1) on_flags (in=in2);
  by event_id;
  if in1 and in2 then output splice_info_w_flags2;
  else if in1 then do;
    flag_cd4_on=0;
    flag_Cd8_on=0;
    flag_cd19_on=0;
    output splice_info_w_flags2;
  end;
  else if in2 then output oops;
run;


/* Add in mean depth */

data depth;
  set splicing.splicing_means_by_celltype;
run;

proc sort data=depth;
   by event_id;
proc sort data=splice_info_w_flags2;
   by event_id;
run;

data splice_info_w_means oops;
  merge splice_info_w_flags2 (in=in1) depth (in=in2);
  by event_id;
  if in1 then output splice_info_w_means;
  else if in2 then output oops;
run;

*if event is "off" then set counts to "n.d.";
data splice_info_w_means2;
  set splice_info_w_means;
  if flag_cd4_on=0 then mean_depth_CD4="n.d.";
  if flag_cd8_on=0 then mean_depth_CD8="n.d.";
  if flag_cd19_on=0 then mean_depth_CD19="n.d.";
run;

/* Get constrast P and FDR values */

data list_of_events;
  set splicing.splicing_results_clean;
  where flag_anova_fdr_05 ne .;
run;

data results;
   set splicing.splicing_results_w_annot_fdr;
   keep event_id CD4_CD8_P CD4_CD19_P CD8_CD19_P ProbF
                 CD4_CD8_FDR CD4_CD19_FDR CD8_CD19_FDR anova_fdr_p;
   rename ProbF=anova_P;
   run;

proc sort data=list_of_events;
  by event_id;
proc sort data=results;
  by event_id;
run;

data results2;
  merge list_of_events (in=in1) results (in=in2);
   by event_id;
  if in1 and in2;
run;


/* Get LS means */

data lsmeans;
  set splicing.anova_lsmeans_diffs;
run;

proc sort data=lsmeans;
  by event_id;
run;

data lsmeans2;
   merge list_of_events (in=in1) lsmeans (in=in2);
   by event_id;
   if in1 and in2;
run;

/* Merge and make summary table */

proc sort data=results2;
   by event_id;
proc sort data=lsmeans2;
   by event_id;
proc sort data=splice_info_w_means2;
   by event_id;
run;

data splice_info_w_results;
   merge splice_info_w_means2 (in=in1) results2 (in=in2);
   by event_id;
   if in1 and in2 then output;
   else if in1 then do; *just to make sure!;
      flag_cd4_on=0;
      flag_cd8_on=0;
      flag_cd19_on=0;
      mean_depth_CD4="n.d.";
      mean_depth_CD8="n.d.";
      mean_depth_CD19="n.d.";
      output;
      end;
run;

data results_by_splicing_for_supp;
   merge splice_info_w_results (in=in1) lsmeans2 (in=in2);
   by event_id;
run;


/* Make permenant */

proc sort data=results_by_splicing_for_supp;
   by chr feature1_start feature1_stop feature2_start feature2_stop;
run;

data splicing.results_by_splicing_for_supp;
   retain event_id gene_id flag_autoimmune_gene flag_diabetes_gene
          num_transcripts transcript_id
          chr feature1_id feature1_type feature1_start feature1_stop
          feature2_id feature2_type feature2_start feature2_stop strand
          flag_junction_annotated flag_intron_retention
          flag_exonskip flag_alt_donor flag_alt_acceptor
          flag_cd4_on flag_cd4_on flag_cd4_on
          mean_depth_cd4 mean_depth_cd8 mean_depth_cd19
          lsmeans_CD4 lsmeans_CD8 lsmeans_CD19
          anova_P anova_fdr_p
          magnitude_cd4cd8 CD4_CD8_P CD4_CD8_FDR
          magnitude_cd4cd19 CD4_CD19_P CD4_CD19_FDR
          magnitude_cd8cd19 CD8_CD19_P CD8_CD19_FDR;

   format CD4_CD8_P best32.;
   format CD4_CD19_P best32.;
   format CD8_CD19_P best32.;
   set results_by_splicing_for_supp;
   drop flag_anova_fdr_05;
   rename CD4_CD8_FDR=CD4_CD8_FDR_P
          CD4_CD19_FDR=CD4_CD19_FDR_P
          CD8_CD19_FDR=CD8_CD19_FDR_P;
run;

proc export data=splicing.results_by_splicing_for_supp
     outfile='!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/results_by_splicing_all.csv' dbms=csv replace;
run;
