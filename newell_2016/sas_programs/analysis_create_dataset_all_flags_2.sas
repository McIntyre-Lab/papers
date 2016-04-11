/* Make final dataset with all flags */

libname dmel "!MCLAB/useful_dmel_data/flybase551/sasdata";
libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";


/* this dataset is with apn>0 data */

proc sort data=dmel.fb551_si_fusions_unique_flagged;
  by fusion_id;
  run;

proc sort data=ribo.output_anova_apn0_2;
  by fusion_id;
  run;

proc sort data=ribo.logrpkm_means;
  by fusion_id;
  run;

proc sort data=ribo.output_flag_resids_apn0_2;
  by fusion_id;
  run;

proc sort data=ribo.flag_fdr_contrast_by_fusion_2;
  by fusion_id;
  run;

proc sort data= ribo.on_calls_gt_apn0;
  by fusion_id;
  run;

proc sort data=ribo.output_covtest_apn0_2;
  by fusion_id;
  run;



data ribo.results_by_fusion_new_model;
  merge dmel.fb551_si_fusions_unique_flagged (in=in1) 
	ribo.logrpkm_means
	ribo.output_flag_resids_apn0_2 
    ribo.output_covtest_apn0_2
	ribo.flag_fdr_contrast_by_fusion_2 
	ribo.on_calls_gt_apn0
	;
  by fusion_id;
  if in1;
  run;

/*
proc export data=ribo.results_by_fusion 
	outfile='!MCLAB/arbeitman/arbeitman_ribotag/results_by_fusion2.csv' 
	label dbms=csv replace;
	putnames=yes;
	run;
    */
