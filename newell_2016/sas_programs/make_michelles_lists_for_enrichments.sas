libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";
libname dmel "!MCLAB/useful_dmel_data/flybase551/sas_data";


/* Need to recreate the following lists in order to perform enrichment tests:
    575genesintersectIPandDalton
    1253genesUnionFruABCDaltonandFear
    3121genesmaleIPvsFemaleIPFDR20

The UnionFruABCDaltonandFear doesnt not need to be remade, those lists haven't changed
The 3121 list is the contrast1 list, if the ipmale vs ipfemale is sig different
    This flag is by exon, need it by gene. Will use the new flag 'gene_sexbias_ip_only' when the male is sig.
    Need to make this flag binary
The 575 list is the overlap between these previous two lists. This will be redone.


*/

*First need the results file by gene, and need FBgn;
data gene_results_no_multi_with_exons;
  set ribo.gene_results_no_multi_with_exons;
  keep symbol_cat
  genes_per_fusion chrom exon_gene_id_cat fbgn_cat gene_sex_bias_input gene_compare_enrichment_bysex gene_sex_bias_ip gene_sexbias_ip_only enrich_in_both_sex_diff_in_ip flag_downstream_fru;
  run;

proc sort data=gene_results_no_multi_with_exons nodupkey;
  by symbol_cat;
  run;

proc compare base=ribo.gene_results_no_multi compare=gene_results_no_multi_with_exons ;
run;
*There are 7889 genes here because this includes genes that are not significant, ~3K are signif;

proc freq data=gene_results_no_multi_with_exons;
  tables gene_sexbias_ip_only;
  run;

data genes_results;
  set gene_results_no_multi_with_exons;
  if gene_sexbias_ip_only = 'male'  then flag_ip_male_bias = 1;
    else if gene_sexbias_ip_only = 'female' then flag_ip_male_bias=0;
    else if gene_sexbias_ip_only = 'ns' or 'mixed' then flag_ip_male_bias=-1;
  run;
*16 missing, those are the 16 mixed;
proc freq data=genes_results;
  tables flag_ip_male_bias;
  run;

/*
Flag_ip_male_bias       freq
-1 (notsig)             6,332
0 (female)              331
1 (male)                1,210
*/

data maleIP_biased_1210;
  set genes_results;
  where flag_ip_male_bias=1;
  run;

data femaleIP_biased_331;
  set genes_results;
  where flag_ip_male_bias=0;
  flag_ip_female_bias=1;
  run;


data ribo.genes_maleIP_biased_1210;
  set maleIP_biased_1210;
  rename fbgn_cat=primary_fbgn;
  run;

data ribo.genes_femaleIP_biased_331;
  set femaleIP_biased_331;
  rename fbgn_cat=primary_fbgn;
  run;



/* Need to do the overlap now between the FruABC list and the MaleIPBiased lists*/
proc sort data=ribo.genes_maleIP_biased_1210;
  by primary_fbgn;
  run;
proc sort data=ribo.unionfruabc_daltonandfear_1253;
  by primary_fbgn;
  run;

data overlap;
  merge ribo.genes_maleIP_biased_1210 (in=in1) ribo.unionfruabc_daltonandfear_1253 (in=in2);
  by primary_fbgn;
  if in1 and in2;
  run;
  /* 250 genes overlap these two lists now*/

data ribo.intersectIPandDalton_250;
  set overlap;
  run;

/* At this point, I merge in the motif flags and counts. these were done
 * previously and do not need to be redone 
 The dataset ribo.motif_flags_and_cnts was imported and set in 'motif_analysis_import_v2.sas'*/

*250 overlap genes;
data list250_fruABC;
  merge ribo.intersectIPandDalton_250 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list250_fruabc;
  set list250_fruabc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;

*1210 IPmale biased genes;
data list1210_fruabc;
  merge ribo.genes_maleIP_biased_1210 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1 ;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list1210_fruabc;
  set list1210_fruabc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;

*331 female biased ip genes;
proc sort data=ribo.genes_femaleIP_biased_331;
  by primary_fbgn;
  run;

data list331_fruabc;
  merge ribo.genes_femaleIP_biased_331 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1 ;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list331_fruabc;
  set list331_fruabc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;


* the 1253 dataset shouldnt be any different. ;

* Make datasets permanent ; 
data ribo.intersectIPandDalton_250_motifs;
  set list250_fruabc;
  run;
data ribo.genes_maleIP_biased_1210_motifs;
  set list1210_fruabc;
  run;
data ribo.genes_femaleIP_biased_331_motifs;
  set list331_fruabc;
  run;

