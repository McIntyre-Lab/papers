ods listing; ods html close;
libname splicing '/mnt/data/splicing';
libname splictmp '/mnt/store/diabetes_sandbox/splicing';

/* Combine LS means output */

data splicing.anova_lsmeans;
  set splictmp.splicing_anova_lsmean_: ;
run;

/* Process LS means */

data lsmeans;
  set splicing.anova_lsmeans;
  keep event_id cell_type estimate;
run;

proc sort data=lsmeans;
  by event_id cell_type;
proc transpose data=lsmeans out=lsmeans_sbys;
  by event_id;
  var estimate;
  id cell_type;
run;

data splicing.anova_lsmeans_diffs;
  set lsmeans_sbys;
  if CD4 > CD8 then magnitude_cd4cd8=-CD4/CD8; else magnitude_cd4cd8=CD8/CD4;
  if CD4 > CD19 then magnitude_cd4cd19=-CD4/CD19; else magnitude_cd4cd19=CD19/CD4;
  if CD8 > CD19 then magnitude_cd8cd19=-CD8/CD19; else magnitude_cd8cd19=CD19/CD8;

  keep event_id magnitude_cd4cd8 magnitude_cd4cd19 magnitude_cd8cd19 CD4 CD8 CD19;
  rename CD4=lsmeans_CD4 CD8=lsmeans_CD8 CD19=lsmeans_CD19;
run;


/* Combine means */

data means_by_event;
   set splictmp.splicing_means_by_celltype_: ;
run;


data splicing.splicing_means_by_celltype;
  set means_by_event;
  length mean_str_cd4 $10.;
  length mean_str_cd8 $10.;
  length mean_str_cd19 $10.;
  length sd_str_cd4 $10.;
  length sd_str_cd8 $10.;
  length sd_str_cd19 $10.;
  length mean_depth_CD4_2 $30.;
  length mean_depth_CD8_2 $30.;
  length mean_depth_CD19_2 $30.;

  mean_str_cd4=strip(put(mean_depth_cd4, 10.3));
  sd_str_cd4=strip(put(sd_depth_cd4, 10.3));
  mean_depth_CD4_2=cat(strip(mean_str_cd4), " ± ",strip(sd_str_cd4));

  mean_str_cd8=strip(put(mean_depth_cd8, 10.3));
  sd_str_cd8=strip(put(sd_depth_cd8, 10.3));
  mean_depth_CD8_2=cat(strip(mean_str_cd8), " ± ",strip(sd_str_cd8));

  mean_str_cd19=strip(put(mean_depth_cd19, 10.3));
  sd_str_cd19=strip(put(sd_depth_cd19, 10.3));
  mean_depth_CD19_2=cat(strip(mean_str_cd19), " ± ",strip(sd_str_cd19));

  keep event_id mean_depth_CD4_2 mean_depth_CD8_2 mean_depth_CD19_2;
  rename mean_depth_CD4_2=mean_depth_CD4
         mean_depth_CD8_2=mean_depth_CD8
         mean_depth_CD19_2=mean_depth_CD19 ;
run;




