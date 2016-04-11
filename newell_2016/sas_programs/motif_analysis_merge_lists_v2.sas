/*Now that I have imported/created all the gene lists and motif flags, merge
 * files together and create necessary flags */

libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';

data analysed_genes;
  set genes_results;
  if chrom = 'X' then flag_x =1; else flag_x =0;
  rename fbgn_cat=primary_fbgn;
  run;

proc sort data=analysed_genes;
  by primary_fbgn;
  run;

/* Merge gene lists with the full gene file, create flag for male_biased,
 * fru_regulation, and overlap */

 data list_1253;
   set ribo.unionfruabc_daltonandfear_1253;
     flag_fru_reg=1;    
  drop gene_secondary_identifier organism fru_a_motif_pos fru_a_motif_eval flag_fru_a_multi_motif 
    fru_b_motif_pos fru_b_motif_eval flag_fru_b_multi_motif 
    fru_c_motif_pos fru_c_motif_eval flag_fru_c_multi_motif
    fru_a_motif_cnt fru_b_motif_cnt fru_c_motif_cnt;
  run;

data merge_1253;
  merge analysed_genes (in=in1) list_1253 (in=in2);
  by primary_fbgn;
  if in1 ;
  if flag_fru_reg ne 1 then flag_fru_reg=0;
  run;

data merge_1253;
  set merge_1253;
  array change _numeric_;
    do over change;
    if change=. then change=0;
    end;
run;


/*Merge in the male biased list, add flag for male bias */

proc sort data=ribo.genes_maleip_biased_1210_motifs;
  by primary_fbgn;
  run;

data list_1210;
  set ribo.genes_maleip_biased_1210_motifs;
  flag_male_bias=1;
  drop  fru_a_motif_pos fru_a_motif_eval flag_fru_a_multi_motif 
    fru_b_motif_pos fru_b_motif_eval flag_fru_b_multi_motif 
    fru_c_motif_pos fru_c_motif_eval flag_fru_c_multi_motif
    fru_a_motif_cnt fru_b_motif_cnt fru_c_motif_cnt;
  run;


data merge_1253_1210;
  merge merge_1253 (in=in1) list_1210 (in=in2);
  by primary_fbgn;
  if flag_fru_reg ne 1 then flag_fru_reg = 0;
  if flag_male_bias ne 1 then flag_male_bias = 0;
  run;

proc freq data=merge_1253_1210;
  tables flag_fru_reg;
  run;*1048, not right...?;

proc freq data=merge_1253_1210;
  tables flag_male_bias;
  run; *1210;


/* If flag_fru_reg = 1 and flag_male_bias= 1 then flag_overlap = 1; there should
 * be 250 genes that overlap */

 data merged_all;
   set merge_1253_1210;
   if flag_fru_reg=1 and flag_male_bias=1 then flag_overlap=1;
   else flag_overlap=0;
   run;

proc freq data=merged_all;
  tables flag_overlap;
  run;
*250 overlap!;

data merged_all;
  set merged_all;
  array change _numeric_;
    do over change;
    if change=. then change=0;
    end;
run;


* Merge in the female list;
proc sort data=ribo.genes_femaleip_biased_331_motifs;
  by primary_fbgn;
  run;

data list_331;
  set ribo.genes_femaleip_biased_331_motifs;
  flag_female_bias=1;
  drop  fru_a_motif_pos fru_a_motif_eval flag_fru_a_multi_motif 
    fru_b_motif_pos fru_b_motif_eval flag_fru_b_multi_motif 
    fru_c_motif_pos fru_c_motif_eval flag_fru_c_multi_motif
    fru_a_motif_cnt fru_b_motif_cnt fru_c_motif_cnt;
  run;


data merged_all2;
  merge merged_all (in=in1) list_331 (in=in2);
  by primary_fbgn;
  if flag_female_bias ne 1 then flag_female_bias = 0;
    run;

/* Make dataset permanent for now */

data ribo.motif_enrichment_revision;
  set merged_all2;
  run;
