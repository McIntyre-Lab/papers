/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Need to ID instances were a NIC junction suggests that two mutually-exclusive exons are joined.

   (1) ID constitutive and non-constitutive exons (or better: ID exons with constitutive/non-con splice sites)
   (2) For non-constitutive, identify its flanking constitutive exons (use annotated junctions to help)
   (3) Pairwise comparison of non-constitutive exons by gene, if they share the same flanking
       constitutive exons but none of the same transcripts, then these are a mutually exclusive pair
   (4) Find NIC junctions that link these two exons
   (5) Count instances where the NIC junction is detected/supported (APN 0,2,5)
*/

/*** (1) Identify constitutive and non-constitutive splice sites ***/

/* (a) count transcripts per gene */

data xs2gene;
  set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id transcript_id2;
run;

proc sort data=xs2gene nodup;
   by gene_id transcript_id2;
proc freq data=xs2gene noprint;
   tables gene_id / out=xs_per_gene;
run;

data xs_per_gene2;
  set xs_per_gene;
  xs_per_gene=COUNT;
  drop COUNT PERCENT;
run;

/* (b) count transcripts per donor and acceptor */

data donor2xs;
   set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id stop transcript_id2;
  rename stop=donor_site;
run;

proc sort data=donor2xs nodup;
   by gene_id donor_site  transcript_id2;
proc freq data=donor2xs noprint;
   by gene_id;
   tables donor_site / out=xs_per_donor;
run;

data xs_per_donor2;
  set xs_per_donor;
  xs_per_donor=COUNT;
  drop COUNT PERCENT;
run;

data acceptor2xs;
   set mm10.mm10_exons_w_info;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
  keep gene_id start transcript_id2;
  rename start=acceptor_site;
run;

proc sort data=acceptor2xs nodup;
   by gene_id acceptor_site transcript_id2;
proc freq data=acceptor2xs noprint;
   by gene_id;
   tables acceptor_site / out=xs_per_acceptor;
run;

data xs_per_acceptor2;
  set xs_per_acceptor;
  xs_per_acceptor=COUNT;
  drop COUNT PERCENT;
run;

/* (c) Flag donors and acceptors as constitutive or not */

proc sort data=xs_per_gene2;
   by gene_id;
proc sort data=xs_per_donor2;
   by gene_id;
proc sort data=xs_per_acceptor2;
   by gene_id;
run;

data flag_constit_donor;
  merge xs_per_donor2 (in=in1) xs_per_gene2 (in=in2);
  by gene_id;
  if in1 and in2;
  if xs_per_donor=xs_per_gene then flag_constit_donor=1;
  else flag_constit_donor=0;
run;

data flag_constit_acceptor;
  merge xs_per_acceptor2 (in=in1) xs_per_gene2 (in=in2);
  by gene_id;
  if in1 and in2;
  if xs_per_acceptor=xs_per_gene then flag_constit_acceptor=1;
  else flag_constit_acceptor=0;
run;

/* (d) Identify exons with non-constitutive donors and acceptors */

data alt_donor;
  set flag_constit_donor;
  where flag_constit_donor=0;
  keep gene_id donor_site;
run;

data exons;
  set mm10.mm10_exons_w_info;
  keep gene_id exon_id start stop ;
  rename stop=donor_site start=acceptor_site;
run;

proc sort data=exons;
   by gene_id donor_site;
proc sort data=alt_donor;
   by gene_id donor_site;
run;

data exons_alt_donor;
   merge exons (in=in1) alt_donor (in=in2);
   by gene_id donor_site;
   if in1 and in2;
run;

data alt_acceptor;
  set flag_constit_acceptor;
  where flag_constit_acceptor=0;
  keep gene_id acceptor_site;
run;

proc sort data=exons_alt_donor;
   by gene_id acceptor_site;
proc sort data=alt_acceptor;
   by gene_id acceptor_site;
run;

data exons_alt_acceptor;
   merge exons_alt_donor (in=in1) alt_acceptor (in=in2);
   by gene_id acceptor_site;
   if in1 and in2;
run;  *these are my "alternative" exons;

/* (2) For non-constitutive (alternative) exons, identify its flanking constitutive exons
   using junctions to aid in this */

data junc_sites;
  set evspl.splicing_events_annot_refseq;
  where flag_junction_annotated=1 or num_transcripts > 0;
  length junction_id $50. ;
  junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
  keep junction_id gene_id chr feature1_stop feature2_start;
  rename feature1_stop=donor_site feature2_start=acceptor_site;
run;

/* (a) Constitutive donors */
data constit_donor;
  set flag_constit_donor;
  where flag_constit_donor=1;
  keep gene_id donor_site;
run;

proc sort data=constit_donor nodup;
  by gene_id donor_site;
proc sort data=junc_sites nodup;
  by gene_id donor_site acceptor_site;
run;

data constit_junc_donor;
  merge constit_donor (in=in1) junc_sites (in=in2);
  by gene_id donor_site;
  if in1 and in2;
run;

* Merge in junctions to constitutive donor exons ;

proc sort data=constit_junc_donor;
  by gene_id acceptor_site;
proc sort data=exons_alt_acceptor;
  by gene_id acceptor_site;
run;

data alt_exon2constit_donor;
  merge exons_alt_acceptor (in=in1) constit_junc_donor (in=in2);
  by gene_id acceptor_site;
  if in1 and in2;
  rename junction_id=junc_from_donor_id chr=donor_chr;
run;

/* (b) Constitutive acceptors */
data constit_acceptor;
  set flag_constit_acceptor;
  where flag_constit_acceptor=1;
  keep gene_id acceptor_site;
run;

proc sort data=constit_acceptor nodup;
  by gene_id acceptor_site;
proc sort data=junc_sites nodup;
  by gene_id acceptor_site donor_site;
run;

data constit_junc_acceptor;
  merge constit_acceptor (in=in1) junc_sites (in=in2);
  by gene_id acceptor_site;
  if in1 and in2;
run;

* Merge in junctions to constitutive donor exons ;

proc sort data=constit_junc_acceptor;
  by gene_id donor_site;
proc sort data=alt_exon2constit_donor;
  by gene_id donor_site;
run;

data alt_exon2constit_acceptor;
  merge alt_exon2constit_donor (in=in1) constit_junc_acceptor (in=in2);
  by gene_id donor_site;
  if in1 and in2;
  rename junction_id=junc_from_acceptor_id chr=acceptor_chr;
run; *7434 possible sites;


