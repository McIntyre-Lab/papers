/* For exon fragments detected, how many multi-gene fragments are detected? */

libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';


/* From the set of gene pairs, I want to pick pairs where both genes have a lot of fragments detected
   so that when I examine the data I can more easily decide which transcripts are detected
   Therefore, for each gene I want to calculate the number of fragments detected as a proportion of the total

   I am going to keep fragments at least 27bp in length so I have the set of most informative fragments */


/* Subset multigene fragments */

proc sort data=event.pairs_to_keep nodup;
  by gene_id;
proc sort data=event.multigene_frags_on;
  by gene_id;
run;

data frags_to_keep;
  merge multigene_frags_on (in=in1) pairs_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* merge fragment2gene */

data frag2gene;
  set mm10.mm10_fragment2exon2gene;
  if end-start > 27;
  keep fragment_id gene_id;
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

data event.subset_gene_pairs;
  set subset_genes;
run;

