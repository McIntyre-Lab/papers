/* For exon fragments detected, how many multi-gene fragments are detected? */

libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

* Do for NSCs;

data frags_on;
  set event.flag_fragment_on;
  where flag_fragment_nsc_on=1;
  keep fragment_id;
run;

* I am going to pull only the fragments assigned to 2 genes for simplicity;
data multigene_frags;
  set mm10.mm10_exon_fragment_flagged;
  where num_genes_per_fragment=2;
  keep fragment_id exon_id transcript_id gene_id flag_unique flag_common flag_constitutive;
run;

proc sort data=frags_on;
  by fragment_id;
proc sort data=multigene_frags;
  by fragment_id;
run;

data multigene_frags_on;
  merge frags_on (in=in1) multigene_frags (in=in2);
  by fragment_id;
  if in1 and in2;
run;
*6958 2-gene fragments detected;


/********************************************************************************/

/* I want to look at only gene pairs where each gene is ONLY in that pair

   e.g. If geneID "5111" overlaps with 2 other genes ("5111|5333" and "5111|933") then I don't want to use "5111" */

data gene_list;
  set multigene_frags_on;
  keep gene_id;
run;

proc sort data=gene_list nodup;
  by gene_id;
run;


data gene_list_stack;
  length unstack_gene_id $15.;
  set gene_list;
  do i=1 by 1 while(scan(gene_id,i,"|") ^= '');
     unstack_gene_id=scan(gene_id,i,"|");
     output;
     end;
run;

proc freq data=gene_list_stack noprint;
  tables unstack_gene_id / out=count_pairs;
run;

data genes_to_drop;
  set count_pairs;
  where count > 1;
  keep unstack_gene_id;
run;

proc sort data=genes_to_drop;
  by unstack_gene_id;
proc sort data=gene_list_stack;
  by unstack_gene_id;
run;

data pairs_to_keep;
  merge gene_list_stack (in=in1) genes_to_drop (in=in2);
  by unstack_gene_id;
  if in2 then delete;
  keep gene_id;
run;

/* Subset multigene fragments */

proc sort data=pairs_to_keep nodup;
  by gene_id;
proc sort data=multigene_frags_on;
  by gene_id;
run;

data frags_to_keep;
  merge multigene_frags_on (in=in1) pairs_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* From the set of gene pairs, I want to pick pairs where both genes have a lot of fragments detected
   so that when I examine the data I can more easily decide which transcripts are detected
   Therefore, for each gene I want to calculate the number of fragments detected as a proportion of the total

   I am going to keep fragments at least 27bp in length so I have the set of most informative fragments */

data frag2gene;
  set mm10.mm10_fragment2exon2gene;
  if end-start > 27;
  keep fragment_id gene_id;
run;


data frags_on;
  set event.flag_fragment_on;
  where flag_fragment_nsc_on=1;
  keep fragment_id;
run;


proc sort data=frag2gene;
  by fragment_id;
proc sort data=frags_on;
  by fragment_id;
run;

data frags_on2gene;
  merge frag2gene (in=in1) frags_on (in=in2);
  by fragment_id;
  if in1 and in2 then flag_fragment_on=1;
  else if in1 then flag_fragment_on=0;
  flag_fragment=1; *set this here so I can count fragments;
run;

/* Calc number of fragments, and fragments on */

proc sort data=frags_on2gene;
  by gene_id;
proc means data=frags_on2gene noprint;
  by gene_id;
  var flag_fragment_on flag_fragment;
   output out=frags_per_gene sum(flag_fragment_on)=num_dtct_fragments sum(flag_fragment)=num_fragments;
run;

*calc proportion on;
data frags_per_gene2;
  set frags_per_gene;
  perc_fragments_on=num_dtct_fragments/num_fragments;
run;

*subset genes in gene pairs;

data genes2keep;
  length unstack_gene_id $15.;
  set frags_to_keep;
  do i=1 by 1 while(scan(gene_id,i,"|") ^= '');
     unstack_gene_id=scan(gene_id,i,"|");
     output;
     end;
  keep gene_id unstack_gene_id;
  rename gene_id=gene_id_cat unstack_gene_id=gene_id;
run;

proc sort data=genes2keep nodup;
  by gene_id gene_id_cat;
proc sort data=frags_per_gene2;
  by gene_id;
run;

data frags_per_gene3;
  merge frags_per_gene2 (in=in1) genes2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=frags_per_gene3;
   by descending perc_fragments_on;
run;

*subset if there are more than 10 fragments detected, and more than 80% of fragments are detected;
data subset_genes;
  set frags_per_gene3;
  if num_dtct_fragments > 10 and perc_fragments_on ge 0.8;
run;

proc sort data=subset_genes;
  by gene_id_cat;
run;



/* Looking at the data, there are several interesting pairs, so I am going to subset some of the interesting ones */ 

data first_in_pair second_in_pair;
   set subset_genes;
   by gene_id_cat;
   if first.gene_id_cat then output;
   else output second_in_pair;
   drop _TYPE_ _FREQ_;
run;

data first2;
  set first_in_pair;
  rename gene_id=gene_id_1
         num_dtct_fragments=num_dtct_frag_1
         num_fragments=num_frag_1
         perc_fragments_on=perc_frag_1;
run;


data second2;
  set second_in_pair;
  rename gene_id=gene_id_2
         num_dtct_fragments=num_dtct_frag_2
         num_fragments=num_frag_2
         perc_fragments_on=perc_frag_2;
run;


proc sort data=first2;
  by gene_id_cat;
proc sort data=second2;
  by gene_id_cat;
run;

data gene_pairs;
  merge first2 (in=in1) second2 (in=in2);
  by gene_id_cat;
  if in1 and in2;
run;

/* Drop pairs where the number of fragments is the same, and/or the gene_id is the same */

data gene_pairs2;
  set gene_pairs;
  if gene_id_1=gene_id_2 then delete;
  if num_frag_2=num_frag_1 and num_dtct_frag_2=num_dtct_frag_1 then delete;
run;

/* Let's take the list of remaining genes and determine how many have at least one unique fragment >27bp in length that is detected */


*stack geneIDs;
data gene1;
  set gene_pairs2;
  keep gene_id_1;
  rename gene_id_1=gene_id;
run;

data gene2;
  set gene_pairs2;
  keep gene_id_2;
  rename gene_id_2=gene_id;
run;

data gene3;
  set gene1 gene2;
run;

data uniq_frag2gene;
  set mm10.mm10_exon_fragment_flagged;
  where flag_unique=1;
  if fragment_end-fragment_start < 27 then delete;
  keep gene_id fragment_id;
run;

proc sort data=frags_on;
  by fragment_id;
proc sort data=uniq_frag2gene;
  by fragment_id;
run;

data uniq_frag_on2gene;
  merge uniq_frag2gene (in=in1) frags_on (in=in2);
  by fragment_id;
  if in1 and in2;
run;

proc sort data=gene3 nodup;
  by gene_id;
proc sort data=uniq_frag_on2gene;
  by gene_id;
run;

data genes_in_pairs_w_uniq_on;
  merge gene3 (in=in1) uniq_frag_on2gene (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id;
run;

proc sort data=genes_in_pairs_w_uniq_on nodup;
  by gene_id;
run; *645 genes;

data gene_uniq1;
  set genes_in_pairs_w_uniq_on;
  rename gene_id=gene_id_1;
run;

data gene_uniq2;
  set genes_in_pairs_w_uniq_on;
  rename gene_id=gene_id_2;
run;

proc sort data=gene_uniq1;
  by gene_id_1;
proc sort data=gene_pairs2;
  by gene_id_1;
run;

data gene_pairs_flag1;
  merge gene_pairs2 (in=in1) gene_uniq1 (in=in2);
  by gene_id_1;
  if in2 then flag_gene_1_uniq=1; else flag_gene_1_uniq=0;
  if in1 then output;
run;


proc sort data=gene_uniq2;
  by gene_id_2;
proc sort data=gene_pairs_flag1;
  by gene_id_2;
run;

data gene_pairs_flag2;
  merge gene_pairs_flag1 (in=in1) gene_uniq2 (in=in2);
  by gene_id_2;
  if in2 then flag_gene_2_uniq=1; else flag_gene_2_uniq=0;
  if in1 then output;
run;

/* Subset genepairs where both genes have a unique piece */

data gene_pairs_subset_uniq;
  set gene_pairs_flag2;
  if flag_gene_1_uniq=1 and flag_gene_2_uniq=1;
run;


