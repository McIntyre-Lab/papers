ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Exclude events (junctions, exon fragments) and transcripts from genes with any multigene components 

First: flag genes that have multi-gene components
*/

/* Get list of genes that have ANY multigene exonic regions  */

data genes_with_multi_fus;
   set mm10.mm10_refseq_fusion_si_info_v2;
   where flag_multigene=1;
   keep primary_gene_id;
   rename primary_gene_id=gene_id;
run;

proc sort data=genes_with_multi_fus nodup;
  by gene_id;
run;

/* Get list of genes that have ANY multigene fragments  */

data genes_with_multi_frag;
   length gene_id2 $15.;
   set mm10.mm10_refseq_exon_fragment_info;
   where num_genes_per_fragment > 1;
   do i=1 by 1 while(scan(gene_id,i,"|") ^= "");
       gene_id2=scan(gene_id,i,"|");
       output;
       end;
   keep gene_id2;
   rename gene_id2=gene_id;
run;

proc sort data=genes_with_multi_frag nodup;
  by gene_id;
run;

/* Get list of genes that have ambiguous (multigene) junctions */

data genes_w_multi_jct;
   length gene_id $15.;
   set evspl.flag_event_redundant_seq;
   where flag_redundant_sequence=1;
   gene_id=scan(event_id,1,":");
   keep gene_id;
run;

proc sort data=genes_w_multi_jct nodup;
   by gene_id;
run;



/* Make list of genes to drop */

data all_genes;
  set mm10.mm10_exon2gene;
  keep gene_id;
run;

proc sort data=all_genes nodup;
  by gene_id;
run;

data genes_w_multigene_seq;
   merge all_genes (in=in1) genes_with_multi_fus (in=in2)
         genes_with_multi_frag (in=in3) genes_w_multi_jct (in=in4);
   by gene_id;
   if in2 then flag_multigene_fusion=1; else flag_multigene_fusion=0;
   if in3 then flag_multigene_fragment=1; else flag_multigene_fragment=0;
   if in4 then flag_ambig_junction=1; else flag_ambig_junction=0;
   if in1 then output;
run;

proc freq data=genes_w_multigene_seq noprint;
   tables flag_multigene_fusion*flag_multigene_fragment*flag_ambig_junction / out=multiflags;
run;

proc print data=multiflags;
run;

/*
        flag_multigene_    flag_multigene_    flag_ambig_
 Obs         fusion            fragment         junction     COUNT    PERCENT

  1            0                  0                0         27352    70.9833
  2            0                  0                1          1115     2.8936
  3            1                  0                0           359     0.9317
  4            1                  0                1             2     0.0052
  5            1                  1                0          7470    19.3860
  6            1                  1                1          2235     5.8002

Okay, check a few of the instances where a gene has a multigene fusion, but no multigene fragments
This is the general scenario:

GENE A:   [=================]
GENE B:                      [================]

i.e. exon of gene B starts 1bp after end of exon of gene A.
Going to use the multigene fusion flag to remove genes, since this can interfere with IR events later

This leaves us with 28467 genes to test, and 10066 to exclude (~26%)
*/

data event.genes_flag_multigene;
   set genes_w_multigene_seq;
   keep gene_id flag_multigene_fusion;
   rename flag_multigene_fusion=flag_multigene_region;
run;

