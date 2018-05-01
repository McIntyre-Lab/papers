/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Pick ~10 or so genes to simulate transcripts from. My criteria are:
   (1) At least 3 transcripts, but fewer than 10
   (2) At least one NIC junction that joins mutually-exclusive alternative exons together
   (3) No multigene!
*/

data xs2gene;
  set event.feature2xs2gene_nomulti;
  keep gene_id transcript_id;
run;

proc sort data=xs2gene nodup;
  by gene_id transcript_id;
proc freq data=xs2gene noprint;
  tables gene_id / out=xs_per_gene;
run;

data genes2sim;
  set xs_per_gene;
  where count > 2 and count < 10;
run;

data nic_genes;
  set event.nic_mxe_with_support;
  where flag_events_nsc_apn_ge5=1 and flag_donor_acceptor_mxe=1;
  keep gene_id;
run;

proc sort data=nic_genes nodup;
  by gene_id;
proc sort data=genes2sim;
  by gene_id;
run;

data genes2sim2;
  merge genes2sim (in=in1) nic_genes (in=in2);
  by gene_id;
  if in1 and in2;
run; *10 genes;

/* Make gene list permenant */

data event.genes_w_nic_junction_10genes;
  set genes2sim2;
  keep gene_id;
run;


/* For these genes:
   (1) Make an exon BED file for extracting exon sequences
   (2) Make a list of exons per transcript, so I can assemble transcript sequences (will do this manually)
   (3) Export any NIC-MXE junctions (so I know what to join)
*/

data genes2sim3;
  set genes2sim2;
  keep gene_id;
run;

data exon2gene;
  set mm10.mm10_exons_w_info;
  keep chrom start stop strand exon_id gene_id;
run;

data exon2xs;
  set mm10.mm10_exon2transcript;
run;


data nic2keep;
  set event.nic_mxe_with_support;
  where flag_events_nsc_apn_ge5=1 and flag_donor_acceptor_mxe=1;
run;

proc sort data=genes2sim3;
  by gene_id;
proc sort data=exon2gene;
  by gene_id;
proc sort data=nic2keep;
  by gene_id;
run;


data exon2gene_keep;
  merge genes2sim3 (in=in1) exon2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

data nic2sim;
  merge genes2sim3 (in=in1) nic2keep (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id chr strand donor_site acceptor_site junction_id donor_transcript_id acceptor_transcript_id;
run;

proc sort data=exon2gene_keep;
  by exon_id;
proc sort data=exon2xs;
  by exon_id;
run;

data exon2xs_keep;
  merge exon2gene_keep (in=in1) exon2xs (in=in2);
  by exon_id;
  if in1 and in2;
run;

/* Export */

data exon_bed;
  retain chrom start stop exon_id score strand;
  length score $1.;
  set exon2gene_keep;
  score=".";
  keep chrom start stop exon_id score strand;
run;

proc sort data=exon2xs_keep;
   by gene_id transcript_id start stop exon_id;
run;

proc export data=exon_bed outfile="!MCLAB/event_analysis/references/exons_10_genes_to_simulate.bed"
    dbms=tab replace;
    putnames=no;
run;

proc export data=exon2xs_keep outfile="!MCLAB/event_analysis/references/exon2transcript_10_genes_to_simulate.tsv"
    dbms=tab replace;
run;

proc export data=nic2sim outfile="!MCLAB/event_analysis/references/nic_junctions_10_genes_to_simulate.tsv"
    dbms=tab replace;
run;

