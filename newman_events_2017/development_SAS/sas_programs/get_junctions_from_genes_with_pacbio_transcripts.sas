
/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* When I compare STAR/catalog junctions against PB, my false negative/positive rates seem "off"
   I am thinking this is because I am including junctions from all possible genes, whereas only
   a subset of genes have PacBio transcripts (in part because PB sort-of reflects what transcripts 
   are well-expressed, and is probably missing low-expressed transcripts.

   So instead, I am going to limit my comparison of junctions to only those that are found within
   genes that have PacBio coverage, as we really can't say anything about
   junctions from genes that aren't in the PB data */

/* So, subset junctions from genes that have PB transcripts. To do this, I am going to iterate over
   all PB genes and extract junctions that fall within the genomic coordinates of that gene */

data pb_genes;
   set evspl.pacbio_exons;
   keep chrom start stop strand gene_id;
run;

proc sort data=pb_genes nodup;
   by chrom strand gene_id start stop;
run;

data pb_gene_start pb_gene_stop;
  set pb_genes;
  by chrom strand gene_id;
  if first.gene_id then output pb_gene_start;
  if last.gene_id then output pb_gene_stop;
run;

data pb_gene_start2;
   set pb_gene_start;
   keep gene_id chrom strand start;
run;

data pb_gene_stop2;
   set pb_gene_stop;
   keep gene_id chrom strand stop;
run;

proc sort data=pb_gene_start2;
   by gene_id chrom strand;
proc sort data=pb_gene_stop2;
   by gene_id chrom strand;
run;

data pb_gene_coord;
  merge pb_gene_start2 (in=in1) pb_gene_stop2 (in=in2);
  by gene_id chrom strand;
  if in1 and in2;
run;

/* iterate over PacBio gene coordinates and keep any junction (by coordinate) that falls within a PB gene */

data all_junc;
  set event.catalog_pacbio_star_junctions;
  keep chr event_id strand donor_stop acceptor_start;
run;

* set up my dataset to append to;
data junc2keep;
   length chr $20.;
   length event_id $450.;
   length strand $1.;
   format donor_stop best32.;
   format acceptor_start best32.;
   if chr = "" then delete;
run;

%let numGenes = 7836; *total number of PB genes;

%macro getJunc();
  /* declare macro variables */
  %local iterGene pbChr pbStart pbStop pbStrand ;

  %let iterGene=1;
   %do %while (&iterGene. <= &numGenes.);

  %put &iterGene.;

   /* Gene coordinates */
   data _null_;
     set pb_gene_coord (firstobs=&iterGene. obs=&iterGene.); *read in 1 record;
     call symputx("pbChr", strip(chrom));
     call symputx("pbStart", strip(start));
     call symputx("pbStop", strip(stop));
     call symputx("pbStrand", strip(strand));
     stop;
   run;

   /* Extract junctions within gene coordinates and append to junc2keep dataset */
 
   data keepers;
     set all_junc;
     where chr="&pbChr." and strand="&pbStrand.";
     if donor_stop ge &pbStart. and acceptor_start le &pbStop. then output;
   run;

   proc append base=junc2keep data=keepers;
   run;

   /*increment gene iterator */
   %let iterGene=%eval(&iterGene.+1);
 
   %end;
%mend;

%getJunc();

data event.all_juncs_in_pb_regions;
   set junc2keep;
run;

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

  
     
