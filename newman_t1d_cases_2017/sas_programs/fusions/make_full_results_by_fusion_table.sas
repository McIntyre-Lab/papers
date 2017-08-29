libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Full exon results table:

chr fusion_start fusion_stop fusion_id gene_id exon_id flag_multigene
flag_immuno_gene flag_ibase_diabetes
flag_cd19_on flag_cd4_on flag_cd8_onflag_all_on
mean_apn_cd4 (w/ SD) mean_apn_cd8 (w/ SD) mean_apn_cd19 (w/ SD)
contrast_cd4cd8_p contrast_cd4cd8_fdr
contrast_cd4cd19_p contrast_cd4cd19_fdr
contrast_cd8cd19_p contrast_cd8cd19_fdr */

/* Get fusion coordinates */

data fus_coord;
  set hg19.unique_info_fusions_si;
  drop fusion_strand num_genes num_exons;
run;

data fus2exon;
  set hg19.hg19_aceview_fusions_si_info;
  keep fusion_id exon_id;
run;

proc sort data=fus2exon nodup;
  by fusion_id exon_id;
proc freq data=fus2exon noprint;
  tables fusion_id / out=exons_per_fusion;
proc sort data=exons_per_fusion;
  by descending count;
proc print data=exons_per_fusion (obs=1);
run; * max 94 exons per fusion;

data cat_exons;
  array exon[94] $ 39.;
  retain exon1-exon94;
  set fus2exon;
  by fusion_id;
  if first.fusion_id then do;
     call missing(of exon1-exon94);
     records = 0;
  end;
  records + 1;
  exon[records]=exon_id;
  if last.fusion_id then output;
run;

*clean up the output file;
data cat_exons2;
  set cat_exons;
  length cat_exon_id $ 1500.;
  cat_exon_id= catx("|", OF exon1-exon94);
  keep records fusion_id cat_exon_id;
  rename records= num_exons_per_region cat_exon_id=exon_id;
  run;

proc sort data=cat_exons2;
  by fusion_id;
proc sort data=fus_coord;
  by fusion_id;
run;

data fus_info;
  merge fus_coord (in=in1) cat_exons2 (in=in2);
  by fusion_id;
  if in1 and in2;
run;


/* Convert gene-level flags to fusion-level flags */

data ai_genes t1d_genes;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai_genes;
   if flag_immunobase_diabetes_gene=1 then output t1d_genes;
   keep gene_id;
run;

data fus2gene;
  set hg19.hg19_aceview_fusions_si_info;
  keep fusion_id gene_id;
run;

proc sort data=fus2gene nodup;
  by gene_id fusion_id;
proc sort data=ai_genes nodup;
  by gene_id;
proc sort data=t1d_genes nodup;
  by gene_id;
run;

data ai_fusions t1d_fusions;
  merge fus2gene (in=in1) ai_genes (in=in2) t1d_genes (in=in3);
  by gene_id;
  if in1 and in2 then output ai_fusions;
  if in1 and in3 then output t1d_fusions;
  keep fusion_id;
run;

proc sort data=fus_info;
   by fusion_id;
proc sort data=ai_fusions nodup;
   by fusion_id;
proc sort data=t1d_fusions nodup;
  by fusion_id;
run;

data fus_info_w_flags;
  merge fus_info (in=in1) ai_fusions (in=in2) t1d_fusions (in=in3);
  by fusion_id;
  if in2 then flag_autoimmune_gene=1; else flag_autoimmune_gene=0;
  if in3 then flag_diabetes_gene=1; else flag_diabetes_gene=0;
  if in1 then output;
run;

/* Add in "on" flags */

data on_flags;
  set con.fusions_on_gt_apn0;
  drop flag_fusion_on0;
run;

proc sort data=on_flags;
  by fusion_id;
run;

data fus_info_w_flags2;
  merge fus_info_w_flags (in=in1) on_flags (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* For each fusion*cell type, calculate the mean normalized APN and SD
   Drop individuals with low coverage! */

data fus_counts;
  set con.fusion_q3_norm_data_all;
  if Name="2009-PC-0144" then delete;
  if Name="2009-PC-0221" then delete;
  if Name="2009-PC-0235" then delete;
  if Name="2009-PC-0236" then delete;
  if Name="2009-PC-0237" then delete;
  q3_q3_apn=(2 ** log_q3_q3_apn) - 1;
  drop log_q3_q3_apn;
run;

data name2cell;
   set con.design_by_subject_new;
   keep name cell_type;
run;

proc sort data=fus_counts;
  by name;
proc sort data=name2cell nodup;
  by name;
run;

data fus_counts2;
  merge name2cell (in=in1) fus_counts (in=in2);
  by name;
  if in1 and in2;
run;

proc sort data=fus_counts2;
  by fusion_id cell_type;
proc means data=fus_counts2 noprint;
  by fusion_id cell_type;
  var q3_q3_apn;
  output out=mean_apn_by_cell mean(q3_q3_apn)=mean_q3_apn
                              stddev(q3_q3_apn)=sd_q3_apn;
run;

data cd4_counts cd8_counts cd19_counts;
  set mean_apn_by_cell;
  if cell_type = "CD4" then output cd4_counts;
  if cell_type = "CD8" then output cd8_counts;
  if cell_type = "CD19" then output cd19_counts;
run;

%macro counts(cell);
data &cell._counts2;
  set &cell._counts;
  length mean_str $10.;
  length sd_str $10.;
  length mean_q3_apn_&cell. $30.;
  mean_str=strip(put(mean_q3_apn, 10.3));
  sd_str=strip(put(sd_q3_apn, 10.3));
  mean_q3_apn_&cell.=cat(strip(mean_str), " ± ",strip(sd_str));
  keep fusion_id mean_q3_apn_&cell.;
run;

proc sort data=&cell._counts2;
   by fusion_id;
run;
%mend;
%counts(cd4);
%counts(cd8);
%counts(cd19);

data all_counts;
  merge cd4_counts2 (in=in1) cd8_counts2 (in=in2) cd19_counts2 (in=in3);
  by fusion_id;
  if in1 and in2 and in3;
run;

proc sort data=all_counts;
  by fusion_id;
proc sort data=fus_info_w_flags2;
  by fusion_id;
run;

data fus_info_w_counts noinfo;
  merge fus_info_w_flags2 (in=in1) all_counts (in=in2);
  by fusion_id;
  if in1 then output fus_info_w_counts;
  else output noinfo;
run;

*if fusion is "off" then set counts to "n.d.";

data fus_info_w_counts2;
  set fus_info_w_counts;
  if flag_CD4_on=0 then mean_q3_apn_cd4="n.d.";
  if flag_CD8_on=0 then mean_q3_apn_cd8="n.d.";
  if flag_CD19_on=0 then mean_q3_apn_cd19="n.d.";
run;

/* Get constrast P and FDR values */

data results;
  set con.results_by_fusion_final;
  keep fusion_id CD4_CD8_Pvalue CD4_CD8_FDR
                 CD4_CD19_Pvalue CD4_CD19_FDR
                 CD8_CD19_Pvalue CD8_CD19_FDR;
run;

proc sort data=results;
  by fusion_id;
proc sort data=fus_info_w_counts2;
  by fusion_id;
run;

data fus_results_w_info;
  merge fus_info_w_counts2 (in=in1) results (in=in2);
  by fusion_id;
  if in1 and in2;
run;


/* Get LS means */

data lsmeans;
  set con.lsmean_fusions;
  keep fusion_id cell_type estimate;
run;

proc sort data=lsmeans;
  by fusion_id cell_type;
proc transpose data=lsmeans out=lsmeans_sbys;
  by fusion_id;
  var estimate;
  id cell_type;
run;

data calc_magnitude;
  set lsmeans_sbys;
  if CD4 > CD8 then magnitude_cd4cd8=-CD4/CD8; else magnitude_cd4cd8=CD8/CD4;
  if CD4 > CD19 then magnitude_cd4cd19=-CD4/CD19; else magnitude_cd4cd19=CD19/CD4;
  if CD8 > CD19 then magnitude_cd8cd19=-CD8/CD19; else magnitude_cd8cd19=CD19/CD8;
  keep fusion_id magnitude_cd4cd8 magnitude_cd4cd19 magnitude_cd8cd19 CD4 CD8 CD19;
  rename CD4=lsmeans_CD4 CD8=lsmeans_CD8 CD19=lsmeans_CD19;
run;

proc sort data=calc_magnitude;
  by fusion_id;
proc sort data=fus_results_w_info;
  by fusion_id;
run;

data results_by_fusion_for_supp;
  merge fus_results_w_info (in=in1) calc_magnitude ;
  by fusion_id;
  if in1 ;
run;

data con.results_by_fusion_for_supp;
   retain chr fusion_start fusion_stop fusion_id flag_multigene
          gene_id num_exons_per_region exon_id flag_autoimmune_gene
          flag_diabetes_gene flag_CD4_on flag_CD8_on flag_CD19_on
          flag_fusion_all_on0 mean_q3_apn_cd4 mean_q3_apn_cd8
          mean_q3_apn_cd19 lsmeans_CD4 lsmeans_CD8 lsmeans_CD19
          magnitude_cd4cd8 CD4_CD8_Pvalue CD4_CD8_FDR
          magnitude_cd4cd19 CD4_CD19_Pvalue CD4_CD19_FDR
          magnitude_cd8cd19 CD8_CD19_Pvalue CD8_CD19_FDR;
   format CD4_CD8_Pvalue best32.;
   format CD4_CD19_Pvalue best32.;
   format CD8_CD19_Pvalue best32.;
   set results_by_fusion_for_supp;
   rename fusion_start=region_start fusion_stop=region_stop
          fusion_id=exonic_region_id flag_fusion_all_on0=flag_all_on
          CD4_CD8_Pvalue=CD4_CD8_P CD4_CD8_FDR=CD4_CD8_FDR_P
          CD4_CD19_Pvalue=CD4_CD19_P CD4_CD19_FDR=CD4_CD19_FDR_P
          CD8_CD19_Pvalue=CD8_CD19_P CD8_CD19_FDR=CD8_CD19_FDR_P;
run;

proc sort data=con.results_by_fusion_for_supp;
   by chr region_start region_stop;
run;

proc export data=con.results_by_fusion_for_supp
     outfile='!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/results_by_exonic_region_all.csv' dbms=csv replace;
run;

