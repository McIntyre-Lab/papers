libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';


/* Looking at the data, there are several interesting pairs, so I am going to subset some of the interesting ones

I want to subset the gene pairs where each gene has at least one detected fragment (>27bp). This will help me to narrow down the list of candidate pairs to examine */ 

data first_in_pair second_in_pair;
   set event.subset_gene_pairs;
   by gene_id_cat;
   if first.gene_id_cat then output;
   else output second_in_pair;
   drop _TYPE_ _FREQ_;
run;

data first2;
  set first_in_pair;
  rename gene_id=gene_id_1 num_dtct_fragments=num_dtct_frag_1
         num_fragments=num_frag_1 perc_fragments_on=perc_frag_1;
run;

data second2;
  set second_in_pair;
  rename gene_id=gene_id_2 num_dtct_fragments=num_dtct_frag_2
         num_fragments=num_frag_2 perc_fragments_on=perc_frag_2;
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

data frags_on;
  set event.flag_fragment_on;
  where flag_fragment_nsc_on=1;
  keep fragment_id;
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
data event.subset_gene_pairs_w_uniq;
  set gene_pairs_flag2;
  if flag_gene_1_uniq=1 and flag_gene_2_uniq=1;
run;

