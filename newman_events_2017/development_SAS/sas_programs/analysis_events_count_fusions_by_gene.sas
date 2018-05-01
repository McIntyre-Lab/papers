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
  output out=num_fus_on_per_gene sum=num_fusion_on;
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
  if not in2 then num_fusion_on=0;
  if num_fusion_on > 0 then flag_gene_expressed=1;
  else flag_gene_expressed=0;
run;

* flag if gene only has multigene fusions;
data gene_mult gene_nomult;
   set mm10.mm10_fusion_si_info_unique;
   length gene_id2 $20.;
   do i=1 by 1 while(scan(gene_id,i,"|") ^="");
      gene_id2=scan(gene_id,i,"|");
      end;
   if flag_multigene=1 then output gene_mult;
   else output gene_nomult;
   keep gene_id2;
   rename gene_id2=gene_id;
run;

proc sort data=gene_mult nodup;
   by gene_id;
proc sort data=gene_nomult nodup;
   by gene_id;
run;

data gene_only_mult;
  merge gene_mult (in=in1) gene_nomult (in=in2);
  by gene_id;
  if in2 then delete;
run;

proc sort data=gene_only_mult;
  by gene_id;
proc sort data=gene_count_fus_on;
  by gene_id;
run;


/* Make permenant */
data event.flag_gene_expressed;
  merge gene_count_fus_on (in=in1) gene_only_mult (in=in2);
  by gene_id;
  if in2 then flag_gene_only_mult=1;
  else flag_gene_only_mult=0;
run;

/* Count genes with at least one fusion detected */
proc freq data=event.flag_gene_expressed;
   tables flag_gene_expressed*flag_gene_only_mult;
run;

/*

 flag_gene_expressed
           flag_gene_only_mult

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  10855 |    813 |  11668
          |  28.17 |   2.11 |  30.28
          |  93.03 |   6.97 |
          |  28.78 | 100.00 |
 ---------+--------+--------+
        1 |  26864 |      0 |  26864
          |  69.72 |   0.00 |  69.72
          | 100.00 |   0.00 |
          |  71.22 |   0.00 |
 ---------+--------+--------+
 Total       37719      813    38532
             97.89     2.11   100.00



26864 of 38532 genes with at least one non-multigene fusion detected
*/


