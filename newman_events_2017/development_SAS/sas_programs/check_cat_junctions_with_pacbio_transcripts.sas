/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* I have used the genomic coordinates of Pacbio genes to decide what junctions to examine for further analysis.
   I can't be sure that the PacBio genes have the same coordinates as RefSeq genes and thus can't be sure that
   I will have all junctions from all possible genes.

   To verify, I am going to grab the genes for junctions that were kept and check that none of the junctions
   from these genes were missed -- use the junction catalog for this */

data junc_coord_kept;
  set event.all_juncs_in_pb_regions;
  keep chr donor_stop acceptor_start strand;
run;

data junc_cat;
   set evspl.splicing_events_annot_refseq;
   keep chr feature1_stop feature2_start strand event_id gene_id;
   rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;


proc sort data=junc_coord_kept nodup;
   by chr strand donor_stop acceptor_start;
proc sort data=junc_cat;
   by chr strand donor_stop acceptor_start;
run;

data junc_kept_genes;
  merge junc_coord_kept (in=in1) junc_cat (in=in2);
  by chr strand donor_stop acceptor_start;
  if in1 and in2;
run;

data junc_kept_genes2;
  set junc_kept_genes;
  keep gene_id;
run;

proc sort data=junc_kept_genes2 nodup;
  by gene_id;
proc sort data=junc_cat;
  by gene_id;
run;

data junc_cat_to_check;
  merge junc_kept_genes2 (in=in1) junc_cat (in=in2);
  by gene_id;
  if in1 and in2;
run;

data junc_cat_to_check2;
  set junc_cat_to_check;
  drop event_id gene_id;
run;


data junc_kept;
  set event.all_juncs_in_pb_regions;
  drop event_id;
run;

proc sort data=junc_cat_to_check nodup;
  by chr strand donor_stop acceptor_start;
proc sort data=junc_kept nodup;
  by chr strand donor_stop acceptor_start;
run;

data junc_kept2;
  merge junc_kept (in=in1) junc_cat_to_check;
  by chr strand donor_stop acceptor_start;
  if in1 then flag_in_pb_coord=1;
  else flaG_in_pb_coord=0;
run;

/* Now merge with big comparison dataset */

data junc_table;
  set event.catalog_pacbio_star_junctions;
run;

data junc_kept3;
  set junc_kept2;
  drop event_id gene_id;
run;

proc sort data=junc_kept3 nodup;
   by chr strand donor_stop acceptor_start;
proc sort data=junc_table;
  by chr strand donor_stop acceptor_start;
run;

data junc_table_pb;
  merge junc_table (in=in1) junc_kept3 (in=in2);
  by chr strand donor_stop acceptor_start;
  if in1 and in2;
run;

/* Make permenant so I can count stuff */

data event.catalog_pacbio_star_junc_pb;
   set junc_table_pb;
run;

