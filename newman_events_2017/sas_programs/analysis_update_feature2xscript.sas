ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Update the assignment of features to genes and transcripts by removing transcripts of genes
   that are not expressed. For now I am going to remove genes that only have multigene fusions
   as these are tricky to handle */

data feat2xs2gene;
  set event.feature2xs2gene;
run;

data xscript2keep;
  set event.flag_xscript_w_gene_on;
  where flag_xscript_gene_exp=1;
  keep transcript_id;
run;

proc sort data=feat2xs2gene;
  by transcript_id;
proc sort data=xscript2keep;
  by transcript_id;
run;

data feat2xs2gene_drop_noexp;
  merge feat2xs2gene (in=in1) xscript2keep (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Make permenant */

data event.feature2xs2gene_exp_only;
  set feat2xs2gene_drop_noexp;
run;

