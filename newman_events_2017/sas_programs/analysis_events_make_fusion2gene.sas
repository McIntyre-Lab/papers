ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Using the set of non-multigene exons, determine what genes are likely expressed.
   Assign exons(fusions)  to gene, removing multigene exons(fusions) */

data fus_info;
   set event.features_w_annotations;
   where feature_type="fusion" and flag_multigene=0;
   keep feature_id flag_feature_on;
run;

data feat2gene;
   set event.feature2xs2gene;
   drop transcript_id;
run;

proc sort data=fus_info;
  by feature_id;
proc sort data=feat2gene nodup;
  by feature_id gene_id;
run;

data single_fus2gene;
  merge fus_info (in=in1) feat2gene (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* Count number of "on" fusions per gene */

proc sort data=single_fus2gene;
  by gene_id;
proc means data=single_fus2gene noprint;
  by gene_id;
  var flag_feature_on;
  output out=num_fus_on_per_gene sum=num_fusions_on;
run;

/* Get list of genes with at least one fusion on */

data gene_list;
   set event.feature2xs2gene;
   keep gene_id;
run;

proc sort data=gene_list nodup;
  by gene_id;
proc sort data=num_fus_on_per_gene (drop=_TYPE_ _FREQ_);
  by gene_id;
run;

data gene_count_fus_on;
  merge gene_list (in=in1) num_fus_on_per_gene (in=in2);
  by gene_id;
  if not in2 then num_fusions_on=0;
  if num_fusions_on > 0 then flag_gene_expressed=1;
  else flag_gene_expressed=0;
run;


/* Make permenant */
data event.flag_gene_expressed;
  set gene_count_fus_on;
run;

/* Count genes with at least one fusion detected */
proc freq data=event.flag_gene_expressed;
   tables flag_gene_expressed;
run;

/*

   flag_gene_                             Cumulative    Cumulative
    expressed    Frequency     Percent     Frequency      Percent
------------------------------------------------------------------
            0       14433       37.46         14433        37.46
            1       24099       62.54         38532       100.00

24099 of 38532 genes with at least one non-multigene fusion detected
*/


