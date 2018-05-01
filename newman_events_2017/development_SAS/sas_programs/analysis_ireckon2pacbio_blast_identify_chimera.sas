/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* iReckon-to-PacBio BLAST: how many chimeric transcripts do we have?

Try the following:
(1) Count transripts with BLAST hits to multiple genes
(2) Count transripts with BLAST hits to multiple genes on the same chromosome
(3) Count transripts with BLAST hits to genes known to overlap (look at fusions)
*/

data blast_hits;
  set event.ireckon2pacbio_blast;
run;

proc sort data=blast_hits;
   by sample query_id transcript_id;
proc freq data=blast_hits noprint;
   by sample query_id;
   tables transcript_id / out=ref_hits_per_query;
run;

data fragmented;
  set ref_hits_per_query;
  where count > 1;
  keep sample query_id transcript_id;
run;

proc sort data=blast_hits;
   by sample query_id transcript_id;
proc sort data=fragmented;
   by sample query_id transcript_id;
run;

data blast_hits_fragmented;
  merge blast_hits (in=in1) fragmented (in=in2);
  by sample query_id transcript_id;
  if in1 and in2;
run;

/* drop hits with more than 5 mismatches, or have gaps */

data drop_gapped_mm_gt5;
  set blast_hits_fragmented;
  if n_mismatch > 5 then delete;
  if n_gapopen > 0 then delete;
run;

/* Get geneID */

data add_geneID;
  length gene_id $15.;
  set drop_gapped_mm_gt5;
  gene_id=catt("PB.",scan(transcript_id,2,"."));
run;

/* Count: number of transcripts that BLAST to multiple transcripts/genes */

data query2gene;
   set add_geneID;
   keep sample query_id gene_id;
run;

data query2xscript;
   set add_geneID;
   keep sample query_id transcript_id;
run;

proc sort data=query2gene nodup;
  by sample query_id gene_id ;
proc freq data=query2gene noprint;
  by sample ;
  tables query_id / out=genes_per_hit;
run;

proc sort data=query2xscript nodup;
  by sample query_id transcript_id ;
proc freq data=query2xscript noprint;
  by sample ;
  tables query_id / out=xscripts_per_hit;
run;

data multigene;
  set genes_per_hit;
  where count > 1;
  keep query_id sample;
run;

data multixscript;
  set xscripts_per_hit;
  where count > 1;
  keep query_id sample;
run;

/* No iReckon transcripts are chimeric RNAs */



