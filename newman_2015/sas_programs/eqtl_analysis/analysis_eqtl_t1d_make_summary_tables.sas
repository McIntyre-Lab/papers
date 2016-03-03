/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname splice '/mnt/data/splice';

/* Need to do:
   How many index SNPs have tested eQTLs (total, exons, junctions, IR)?
       Number per index SNP, number of unique SNPs per index SNP
   How many index SNPs have sig eQTLs (total, exons, junctions, IR)?
       Number per index SNP, number of unique SNPs per index SNP
*/

/* Number of index SNPs with tested eQTLs */

data index_snp;
   set eqtl.eqtl_results_w_onengut_index;
   keep onengut_index_snp;
run;
proc sort data=index_snp nodup;
   by onengut_index_snp;
run;
* 50 index SNPs with SNPs tested;

/* Number of index SNPs with tested eQTLs (R2>0.8) */
data index_snp;
   set eqtl.eqtl_results_w_onengut_index;
   where onengut_snp2eqtl_snp_r2 gt 0.8;
   keep onengut_index_snp;
run;
proc sort data=index_snp nodup;
   by onengut_index_snp;
run;
* 41 index SNPs with SNPs tested and r2 gt 0.8;


/* Number of index SNPs with sig eQTLs */
data index_snp;
   set eqtl.eqtl_results_w_onengut_index;
   if flag_eqtl_sig=1;
   keep onengut_index_snp;
run;
proc sort data=index_snp nodup;
   by onengut_index_snp;
run;
* 37 index SNPs with significant eQTL tested;

/* Number of index SNPs with sig eQTLs (R2>0.8) */
data index_snp;
   set eqtl.eqtl_results_w_onengut_index;
   if flag_eqtl_sig=1 and onengut_snp2eqtl_snp_r2 gt 0.8;
   keep onengut_index_snp;
run;
proc sort data=index_snp nodup;
   by onengut_index_snp;
run;
* 18 index SNPs with SNPs tested and r2 gt 0.8;

/* Make summary tables */

/* Need to get coordinates for features */

data exons_tested splicing_tested;
  set eqtl.eqtl_results_w_onengut_index;
  if feature_type='exon' then output exons_tested;
  else output splicing_tested;
  keep feature_id;
run;

proc sort data=exons_tested nodup;
   by feature_id;
proc sort data=splicing_tested nodup;
   by feature_id;
run;

data exon_coordinates;
   set fus.unique_info_fusions_si;
   length feature_pos $50.;
   feature_pos=cats('chr',chr,':',fusion_start,'-',fusion_stop);
   keep fusion_id feature_pos;
   rename fusion_id=feature_id;
run;


data splicing_coordinates;
   set splice.splicing_events_annotations;
   length feature_pos $50.;
   if feature1_size ge 37 then feature1_start2=feature1_stop-37;
   else  feature1_start2=feature1_start;
   if feature2_size ge 37 then feature2_stop2=feature2_start+37;
   else feature2_stop2=feature2_stop;
   feature2_length2=feature1_stop-feature1_start2;
   feature1_length2=feature2_stop2-feature2_start;
   if flag_intron_retention=1 then do;
      feature_pos=cats('chr',chr,':',feature1_start2,'-',feature2_stop2);
      end;
   else do;
      feature_pos=cats('chr',chr,':',feature1_start2,'-',feature1_stop,':',feature2_start,'-',feature2_stop2);
      end;
   keep event_id feature_pos feature1_length2 feature2_length2; 
   rename event_id=feature_id;
run;

proc sort data=exon_coordinates;
   by feature_id;
proc sort data=splicing_coordinates;
   by feature_id;
run;


data exons_tested_coord;
   merge exons_tested (in=in1) exon_coordinates (in=in2);
   by feature_id;
   if in1 and in2;
run;

data splicing_tested_coord;
   merge splicing_tested (in=in1) splicing_coordinates (in=in2);
   by feature_id;
   if in1 and in2;
run;

data feature_coord;
  set splicing_tested_coord exons_tested_coord;
run;

proc sort data=feature_coord;
   by feature_id;
proc sort data=eqtl.eqtl_results_w_onengut_index;
   by feature_id;
run;

data eqtl.eqtl_results_w_onengut_index2 no_Eqtl_oops no_coord_oops;
   merge eqtl.eqtl_results_w_onengut_index (in=in1) feature_coord (in=in2);
   by feature_id;
   if in1 and in2 then output eqtl.eqtl_results_w_onengut_index2;
   else if in1 then output no_coord_oops;
   else output no_eqtl_oops;
   drop feature2_length2 feature1_length2;
run;

/* First table made. Now onto the second table. What I want to do:
   Per index_snp and tested SNP
   1. bin credible SNPs into different r2 ranges
   2. cat together credible SNP id and r2 value
   3. cat together all credible SNPs within bins
   Essentially, first table with the credible SNPs put together
*/


data eqtl_results_w_og_index_r2_bins;
    set eqtl.eqtl_results_w_onengut_index2;
    length r2 $5.;
    length r2_bin $3.;
    length credible_snp_w_r2 $25.;
    if onengut_snp2eqtl_snp_r2 lt 0.01 then r2 = '<0.01';
    else if onengut_snp2eqtl_snp_r2 eq 1 then r2 = '1.00';
    else do;
         r2_temp=round(onengut_snp2eqtl_snp_r2, 0.01);
         r2=put(r2_temp, 4.2); end;
    credible_snp_w_r2=catt(onengut_snp_id,' (',r2,')');
    /* Add in r2 bins */
    if onengut_snp2eqtl_snp_r2 ge 0.8 then r2_bin="80%";
    else if onengut_snp2eqtl_snp_r2 ge 0.5 then r2_bin="50%";
    else r2_bin="Low";
    drop r2_temp;
run;


* cat together duplicate credible SNPs;
* First need the max number of credible SNPs per bin per eQTL test;

data eqtl_results_w_og_index_r2_bins2;
   set eqtl_results_w_og_index_r2_bins;
   keep snp_id feature_id credible_snp_w_r2 r2_bin;
run;

proc sort data=eqtl_results_w_og_index_r2_bins2 nodup;
   by snp_id feature_id credible_snp_w_r2;
run;

proc freq data=eqtl_results_w_og_index_r2_bins2 noprint;
   where r2_bin="80%";
   by snp_id;
   tables feature_id / out=cred_snp_count_80perc;
run;

proc sort data=cred_snp_count_80perc;
  by descending count;
run;

proc freq data=eqtl_results_w_og_index_r2_bins2 noprint;
   where r2_bin="50%";
   by snp_id;
   tables feature_id / out=cred_snp_count_50perc;
run;

proc sort data=cred_snp_count_50perc;
  by descending count;
run;

proc freq data=eqtl_results_w_og_index_r2_bins2 noprint;
   where r2_bin="Low";
   by snp_id;
   tables feature_id / out=cred_snp_count_low;
run;

proc sort data=cred_snp_count_low;
  by descending count;
run;

*19 credible SNPs max for r2>80%;
*16 credible SNPs max for r2>50%;
*18 credible SNPs max for r2<50%;

%macro cat_r2_snps(cutoff, cutoff_num);
data cat_og_snps_&cutoff_num.; 
  array og_snp[19] $ 25.;
  retain og_snp1-og_snp19;
  set eqtl_results_w_og_index_r2_bins2;
  where r2_bin="&cutoff.";
  by snp_id feature_id;
  if first.feature_id then do;
     call missing(of og_snp1-og_snp19);
     records = 0;
  end;
  records + 1;
  og_snp[records]=credible_snp_w_r2;
  if last.feature_id then output;
run;

  *clean up the output file;

data cat_og_snps_&cutoff_num._2;
  set cat_og_snps_&cutoff_num.; 
  length credible_snp_r2_&cutoff_num. $ 600.;
         credible_snp_r2_&cutoff_num.= catx(", ", OF og_snp1-og_snp19);
  keep snp_id feature_id credible_snp_r2_&cutoff_num. records;
  rename records=num_credible_snps_&cutoff_num.perc;
  run;
%mend;

%cat_r2_snps(80%, 80);
%cat_r2_snps(50%, 50);
%cat_r2_snps(Low, low);

/* Cat together O-G index SNPs */

data eqtl_index_snps_for_cat;
  set eqtl_results_w_og_index_r2_bins; 
  keep snp_id feature_id onengut_index_snp;
run;

* First un-cat existing index snp IDs;

data eqtl_index_snps_for_cat2;
   length onengut_index_snp2 $15.; 
   set eqtl_index_snps_for_cat;
   do i=1 by 1 while(scan(onengut_index_snp,i,'|') ^=' ');
      onengut_index_snp2=scan(onengut_index_snp,i,'|');
      output;
   end;
   drop i onengut_index_snp;
   rename onengut_index_snp2=onengut_index_snp;
run;

proc sort data=eqtl_index_snps_for_cat2 nodup;
   by snp_id feature_id onengut_index_snp;
run;

proc freq data=eqtl_index_snps_for_cat2 noprint;
   by snp_id;
   tables feature_id / out=index_snp_count;
run;

proc sort data=index_snp_count;
   by descending count;
run;

*3 index SNPs max;

data cat_og_index_snps; 
  array index_snp[3] $ 25.;
  retain index_snp1-index_snp3;
  set eqtl_index_snps_for_cat2;
  by snp_id feature_id;
  if first.feature_id then do;
     call missing(of index_snp1-index_snp3);
     records = 0;
  end;
  records + 1;
  index_snp[records]=onengut_index_snp;
  if last.feature_id then output;
run;

  *clean up the output file;

data cat_og_index_snps_2;
  set cat_og_index_snps; 
  length onengut_index_snp_cat $ 120.;
         onengut_index_snp_cat= catx(", ", OF index_snp1-index_snp3);
  keep snp_id feature_id onengut_index_snp_cat;
  run;

/* Merge concatentated SNP lists with eQTL summary data */

data eqtl_summary_data;
  set eqtl_results_w_og_index_r2_bins;
  drop onengut_snp_id onengut_index_snp onengut_snp2eqtl_snp_r2 r2 r2_bin credible_snp_w_r2;
run;

proc sort data=eqtl_summary_data nodup;
   by snp_id feature_id;
proc sort data=cat_og_index_snps_2;
   by snp_id feature_id;
proc sort data=cat_og_snps_80_2;
   by snp_id feature_id;
proc sort data=cat_og_snps_50_2;
   by snp_id feature_id;
proc sort data=cat_og_snps_low_2;
   by snp_id feature_id;
run;

data eqtl_results_w_onengut_credible oops;
   merge eqtl_summary_data (in=in1) cat_og_index_snps_2 (in=in2) cat_og_snps_80_2 (in=in3) cat_og_snps_50_2 (in=in4)  cat_og_snps_low_2 (in=in5) ;
   by snp_id feature_id;
   if in1 and in2 then do;
      if not in3 then num_credible_snps_80perc=0;
      if not in4 then num_credible_snps_50perc=0;
      if not in5 then num_credible_snps_lowperc=0;
       output eqtl_results_w_onengut_credible; end;
   else output oops;
run;

/* Make permenant */

data eqtl.eqtl_results_w_onengut_credible;
  retain snp_id onengut_index_snp_cat gene_id flag_diabetes_gene feature_id feature_type feature_pos flag_multigene
  flag_eqtl_sig CD4_FDR_P CD8_FDR_P CD19_FDR_P
  num_credible_snps_80perc credible_snp_r2_80 num_credible_snps_50perc credible_snp_r2_50
  num_credible_snps_lowperc credible_snp_r2_low;
  set eqtl_results_w_onengut_credible;
run;

/* Third table. What I want is summaries per index SNP. What I want:
O-G Index SNP
O-G Index SNP credible gene?
Number of sig eQTLs where eQTL SNP is r2>0.8 with index SNP
   exons, junctions, ir
Number of sig eQTLs where eQTL SNP is r2>0.5 with index SNP
Number of sig eQTLs where eQTL SNP is r2<0.5 with index SNP
List of genes from eQTLs (r2>0.8 only)
List/Number of credible SNPs linked (r2>0.8 only)
*/

/* Step 1. Get the list of O-G Index SNPs and credible genes */

data og_index2gene;
   set eqtl.onengut_supptable1;
   keep chr index_snp_position index_snp_rs genes;
run;

proc sort data=og_index2gene nodup;
   by index_snp_rs genes;
run;

* cat together genes;
* how many per SNP?;
proc freq data=og_index2gene noprint;
   tables index_snp_rs / out=index_count;
run;

proc sort data=index_count;
   by descending count;
run;
*8 genes max;

data cat_og_index_genes; 
  array gene[8] $ 12.;
  retain gene1-gene8;
  set og_index2gene;
  by index_snp_rs;
  if first.index_snp_rs then do;
     call missing(of gene1-gene4);
     records = 0;
  end;
  records + 1;
  gene[records]=genes;
  if last.index_snp_rs then output;
run;

  *clean up the output file;

data cat_og_index_genes_2;
  set cat_og_index_genes; 
  length og_index_snp_genes $ 120.;
         og_index_snp_genes= catx(", ", OF gene1-gene8);
  keep chr index_snp_position index_snp_rs og_index_snp_genes;
  run;


/* Step 2. Number of sig eQTLs */

data sig_eqtl_r2_80 sig_eqtl_r2_50 sig_eqtl_r2_low; 
   set eqtl.eqtl_results_w_onengut_index2;
   if flag_eqtl_sig=1;
   if onengut_snp2eqtl_snp_r2 ge 0.8 then output sig_eqtl_r2_80;
   else if onengut_snp2eqtl_snp_r2 ge 0.5 then output sig_eqtl_r2_50;
   else output sig_eqtl_r2_low;
   keep gene_id feature_id feature_type snp_id onengut_index_snp;
run;

/* Get counts for the number of sig eQTLs */

proc freq data=sig_eqtl_r2_80 noprint;
   tables onengut_index_snp / out=sig_eqtl_r2_80_all;
run;

proc freq data=sig_eqtl_r2_80 noprint;
   where feature_type='exon';
   tables onengut_index_snp / out=sig_eqtl_r2_80_exon;
run;

proc freq data=sig_eqtl_r2_80 noprint;
   where feature_type='Junc';
   tables onengut_index_snp / out=sig_eqtl_r2_80_junc;
run;

proc freq data=sig_eqtl_r2_80 noprint;
   where feature_type='IR';
   tables onengut_index_snp / out=sig_eqtl_r2_80_ir;
run;

proc freq data=sig_eqtl_r2_50 noprint;
   tables onengut_index_snp / out=sig_eqtl_r2_50_all;
run;

proc freq data=sig_eqtl_r2_low noprint;
   tables onengut_index_snp / out=sig_eqtl_r2_low_all;
run;


/* Step 3. List of genes from sig eQTLs (r2>0.8) */

data sig_eqtl_genes;
   set eqtl.eqtl_results_w_onengut_index2;
   if flag_eqtl_sig=1 and onengut_snp2eqtl_snp_r2 ge 0.8;
   keep gene_id onengut_index_snp;
run;

* Stack genes;

data sig_eqtl_genes_stack;
   length gene_id2 $36.; 
   set sig_eqtl_genes;
   do i=1 by 1 while(scan(gene_id,i,'|') ^=' ');
      gene_id2=scan(gene_id,i,'|');
      output;
   end;
   drop i gene_id;
run;

* drop duplicates;
proc sort data=sig_eqtl_genes_stack nodup;
   by onengut_index_snp gene_id2;
run;

* cat genes per index_snp;
proc freq data=sig_eqtl_genes_stack noprint;
   tables onengut_index_snp / out=index_count;
run;

proc sort data=index_count;
   by descending count;
run;
*4 genes max;

data cat_eqtl_genes; 
  array gene[4] $ 36.;
  retain gene1-gene4;
  set sig_eqtl_genes_stack;
  by onengut_index_snp;
  if first.onengut_index_snp then do;
     call missing(of gene1-gene4);
     records = 0;
  end;
  records + 1;
  gene[records]=gene_id2;
  if last.onengut_index_snp then output;
run;

  *clean up the output file;

data cat_eqtl_genes_2;
  set cat_eqtl_genes; 
  length eqtl_gene_list $ 160.;
         eqtl_gene_list= catx(", ", OF gene1-gene4);
  keep onengut_index_snp eqtl_gene_list;
  run;

/* Step 4. Number of linked credible SNPs r2>0.8 */

data sig_eqtl_credible_snps;
   set eqtl.eqtl_results_w_onengut_index2;
   if flag_eqtl_sig=1 and onengut_snp2eqtl_snp_r2 ge 0.8;
   keep onengut_snp_id onengut_index_snp;
run;

proc sort data=sig_eqtl_credible_snps nodup;
    by onengut_index_snp onengut_snp_id;
run;

proc freq data=sig_eqtl_credible_snps noprint;
   tables onengut_index_snp / out=og_credible_snp_count;
run;

/* Step 5. Merge all together */

data sig_eqtl_r2_80_all_2;
   set sig_eqtl_r2_80_all;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_80_all;
run;

data sig_eqtl_r2_80_exon_2;
   set sig_eqtl_r2_80_exon;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_80_exons;
run;

data sig_eqtl_r2_80_junc_2;
   set sig_eqtl_r2_80_junc;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_80_junc;
run;

data sig_eqtl_r2_80_ir_2;
   set sig_eqtl_r2_80_ir;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_80_ir;
run;

data sig_eqtl_r2_50_all_2;
   set sig_eqtl_r2_50_all;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_50_all;
run;

data sig_eqtl_r2_low_all_2;
   set sig_eqtl_r2_low_all;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_r2_low_all;
run;

data og_credible_snp_count_2;
   set og_credible_snp_count;
   keep onengut_index_snp count;
   rename count=num_sig_eqtl_credible_r2_80;
run;

proc sort data=sig_eqtl_r2_80_all_2;
   by onengut_index_snp;
proc sort data=sig_eqtl_r2_80_exons_2;
   by onengut_index_snp;
proc sort data=sig_eqtl_r2_80_junc_2;
   by onengut_index_snp;
proc sort data=sig_eqtl_r2_80_ir_2;
   by onengut_index_snp;
proc sort data=sig_eqtl_r2_50_all_2;
   by onengut_index_snp;
proc sort data=sig_eqtl_r2_low_all_2;
   by onengut_index_snp;
proc sort data=og_credible_snp_count_2;
   by onengut_index_snp;
proc sort data=cat_eqtl_genes_2;
   by onengut_index_snp;
run;

data og_index_sig_eqtls;
   merge sig_eqtl_r2_80_all_2 (in=in1) sig_eqtl_r2_80_exon_2 (in=in2) sig_eqtl_r2_80_junc_2 (in=in3)
   sig_eqtl_r2_80_ir_2 (in=in4) sig_eqtl_r2_50_all_2 (in=in5) sig_eqtl_r2_low_all_2 (in=in6)
   og_credible_snp_count_2 (in=in7) cat_eqtl_genes_2 (in=in8);
   by onengut_index_snp;
   if not in1 then num_sig_eqtl_r2_80_all=0;
   if not in2 then num_sig_eqtl_r2_80_exons=0;
   if not in3 then num_sig_eqtl_r2_80_junc=0;
   if not in4 then num_sig_eqtl_r2_80_ir=0;
   if not in5 then num_sig_eqtl_r2_50_all=0;
   if not in6 then num_sig_eqtl_r2_low_all=0;
   if not in7 then num_sig_eqtl_credible_r2_80=0;
run;


/* Merge in  O-G index snp genes */

data cat_og_index_genes_3;
  set cat_og_index_genes_2;
  rename index_snp_rs=onengut_index_snp;
run;

proc sort data=cat_og_index_genes_3;
   by onengut_index_snp;
proc sort data=og_index_sig_eqtls;
   by onengut_index_snp;
run;

data eqtl_results_by_onengut_index oops;
   merge cat_og_index_genes_3 (in=in1) og_index_sig_eqtls (in=in2);
   by onengut_index_snp;
   if in1 and in2 then output eqtl_results_by_onengut_index;
   else if in1 then do;
       num_sig_eqtl_r2_80_all=0;
       num_sig_eqtl_r2_80_exons=0;
       num_sig_eqtl_r2_80_junc=0;
       num_sig_eqtl_r2_80_ir=0;
       num_sig_eqtl_r2_50_all=0;
       num_sig_eqtl_r2_low_all=0;
       num_sig_eqtl_credible_r2_80=0;
       output eqtl_results_by_onengut_index; end;
   else output oops;
run;


/* Reorder variables and make permenant */

data eqtl.eqtl_results_by_onengut_index;
  retain chr index_snp_position onengut_index_snp og_index_snp_genes
  num_sig_eqtl_r2_80_all num_sig_eqtl_credible_r2_80 eqtl_gene_list
  num_sig_eqtl_r2_80_exons num_sig_eqtl_r2_80_junc num_sig_eqtl_r2_80_ir
  num_sig_eqtl_r2_50_all num_sig_eqtl_r2_low_all;
  set eqtl_results_by_onengut_index;
run;

data eqtl.eqtl_results_by_onengut_index2;
  set eqtl.eqtl_results_by_onengut_index;
  if num_sig_eqtl_r2_80_all gt 0;
  eqtl_gene_list2=tranwrd(eqtl_gene_list, "and", "/");
  keep chr index_snp_position onengut_index_snp num_sig_eqtl_r2_80_all num_sig_eqtl_credible_r2_80 eqtl_gene_list2;
run;


/* EXPORT */

proc export data=eqtl.eqtl_results_by_onengut_index outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_results_by_onengut_index_snp.tsv' dbms=tab replace;
run;

proc export data=eqtl.eqtl_results_by_onengut_index2 outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_results_by_onengut_index_snp_v2.tsv' dbms=tab replace;
run;

proc export data=eqtl.eqtl_results_w_onengut_credible outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_results_with_onengut_credible_snps.tsv' dbms=tab replace;
run;

proc export data=eqtl.eqtl_results_w_onengut_index2 outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_results_with_onengut_index_snps.tsv' dbms=tab replace;
run;

