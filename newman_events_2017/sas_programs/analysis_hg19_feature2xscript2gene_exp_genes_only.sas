libname event '!MCLAB/event_analysis/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splice '/mnt/data/splice';
libname splicing '/mnt/store/splice';

/* Create index of fragments and junctions to transcript and gene */

data junc2xs;
  set splicing.splicing_events_annotations;
  where transcript_id ne '';
  keep event_id gene_id transcript_id;
run;

data frag2xs;
   set hg19.hg19_exon_fragment_flagged;
   keep fragment_id gene_id transcript_id;
run;

data junc2xs_stack;
   length transcript_id2 $60.;
   set junc2xs;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      output;
      end;
   drop transcript_id  i;
   rename transcript_id2=transcript_id;
run;

data frag2xs_stack;
   length transcript_id2 $60.;
   set frag2xs;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      output;
      end;
   drop transcript_id gene_id i;
   rename transcript_id2=transcript_id;
run;

data frag2xs_stack2;
   length fragment_id $2550.;
   length gene_id $36.;
   set frag2xs_stack;
   gene_id=scan(transcript_id,1,".");
   rename fragment_id=feature_id;
run;

data junc2xs_stack2;
   set junc2xs_stack;
   rename event_id=feature_id;
run;

data feature2xs2gene;
   set junc2xs_stack2 frag2xs_stack2;
run;

/* Remove events of genes not expressed, or have multigene parts */

data gene_exp;
  set event.hg19_flag_gene_expressed;
  if flag_multigene=1 then delete;
  if flag_cd4_gene_on=1 or flag_cd8_gene_on=1 or flag_cd19_gene_on=1 then output;
  keep gene_id;
run; *27993 genes kept;

proc sort data=feature2xs2gene;
  by gene_id;
proc sort data=gene_exp;
  by gene_id;
run;

data feature2xs2gene_exp_only;
   merge feature2xs2gene (in=in1) gene_exp (in=in2);
   by gene_id;
   if in1 and in2;
run;

/* Make permenant */

data event.hg19_feature2xs2gene_exp_only;
   set feature2xs2gene_exp_only;
run;

