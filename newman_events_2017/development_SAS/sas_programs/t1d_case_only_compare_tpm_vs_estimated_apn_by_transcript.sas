ods listing; ods html close;
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Now I want to compare APN vs TPM : first for all transcripts, then the MEI per gene per subject */

data tpm_data;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  log_tpm=log(tpm+1);
run;

data est_apn;
   set eventloc.hg19_rsem_reduced_xs_apn_est;
run;

data design;
  set con.design_by_subject_new;
  if name =  '2009-PC-0221' then delete; *sample 75 cd8;
  if name =  '2009-PC-0144' then delete; *sample 48 cd4;
  if name =  '2009-PC-0236' then delete; *sample 80;
  if name =  '2009-PC-0237' then delete; *sample 80;
  if name =  '2009-PC-0235' then delete; *sample 80;
  keep library name subject_id;
run;

proc sort data=est_apn;
  by name;
proc sort data=design;
  by name;
run;

data est_apn2;
  merge design (in=in1) est_apn (in=in2);
  by name;
  if in1 and in2;
run;

proc sort data=est_apn2;
  by library transcript_id;
proc sort data=tpm_data;
  by library transcript_id;
run;

data tpm_vs_apn;
  merge est_apn2 (in=in1) tpm_data (in=in2);
  by library transcript_id;
  if in1 and in2;
run;

/* get MEIs per sample */

data flag_mei;
  set event.t1d_flag_mei_cell_sbys;
  keep gene_id subject_id transcript_id flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
run;

proc sort data=flag_mei;
  by subject_id gene_id transcript_id;
proc transpose data=flag_mei out=stack_mei;
  by subject_id gene_id transcript_id;
  var flag_mei_CD19 flag_mei_CD4 flag_mei_cd8;
run;

data stack_mei2;
  set stack_mei;
  length cell_type $4.;
  if _NAME_ = "flag_mei_CD19" then cell_type="CD19";
  else if _NAME_ = "flag_mei_CD4" then cell_type="CD4";
  else if _NAME_ = "flag_mei_CD8" then cell_type="CD8";
  else delete;
  drop _NAME_;
  rename COL1=flag_mei;
run;

proc sort data=stack_mei2;
  by subject_id cell_type transcript_id;
proc sort data=tpm_vs_apn;
  by subject_id cell_type transcript_id;
run;

data tpm_vs_apn_flag_mei;
  merge tpm_vs_apn (in=in1) stack_mei2 (in=in2);
  by subject_id cell_type transcript_id;
  if in1 and in2;
run;

/* Calc correlation */

proc corr data=tpm_vs_apn_flag_mei pearson;
  var tpm log_tpm apn log_apn;
run;

/*
          Pearson Correlation Coefficients, N = 5595570
                    Prob > |r| under H0: Rho=0

                    TPM       log_tpm           apn       log_apn

  TPM           1.00000       0.30001       0.37185       0.18224
                               <.0001        <.0001        <.0001

  log_tpm       0.30001       1.00000       0.25130       0.43576
                 <.0001                      <.0001        <.0001

  apn           0.37185       0.25130       1.00000       0.46379
                 <.0001        <.0001                      <.0001

  log_apn       0.18224       0.43576       0.46379       1.00000
                 <.0001        <.0001        <.0001

*/

proc corr data=tpm_vs_apn_flag_mei pearson;
  where flaG_mei=1;
  var tpm log_tpm apn log_apn;
run;

/*
              Pearson Correlation Coefficients, N = 706971
                       Prob > |r| under H0: Rho=0

                       TPM       log_tpm           apn       log_apn

     TPM           1.00000       0.46924       0.73545       0.37733
                                  <.0001        <.0001        <.0001

     log_tpm       0.46924       1.00000       0.47779       0.89720
                    <.0001                      <.0001        <.0001

     apn           0.73545       0.47779       1.00000       0.46966
                    <.0001        <.0001                      <.0001

     log_apn       0.37733       0.89720       0.46966       1.00000
                    <.0001        <.0001        <.0001

*/

proc sort data=tpm_vs_apn_flag_mei;
  by cell_type;
run;

proc corr data=tpm_vs_apn_flag_mei pearson;
  by cell_type;
  var log_tpm log_apn;
run;


/* CD19:
               log_tpm       log_apn

 log_tpm       1.00000       0.43329
                              <.0001

 log_apn       0.43329       1.00000
                <.0001

CD4:
       Prob > |r| under H0: Rho=0

                log_tpm       log_apn

  log_tpm       1.00000       0.43873
                               <.0001

  log_apn       0.43873       1.00000
                 <.0001


CD8:
               log_tpm       log_apn

 log_tpm       1.00000       0.43541
                              <.0001

 log_apn       0.43541       1.00000
                <.0001


*/


proc corr data=tpm_vs_apn_flag_mei pearson;
  by cell_type;
  where flaG_mei=1;
  var log_tpm log_apn;
run;


/*
CD19:
                 log_tpm       log_apn

   log_tpm       1.00000       0.90550
                                <.0001

   log_apn       0.90550       1.00000
                  <.0001


CD4:
               log_tpm       log_apn

 log_tpm       1.00000       0.89887
                              <.0001

 log_apn       0.89887       1.00000
                <.0001



CD8:
                log_tpm       log_apn

  log_tpm       1.00000       0.88708
                               <.0001

  log_apn       0.88708       1.00000
                 <.0001


*/

/* Export counts */

data export_counts;
   set tpm_vs_apn_flag_mei;
   keep subject_id cell_type apn log_apn tpm log_tpm transcript_id flag_mei;
run;

proc export data=export_counts outfile="!MCLAB/event_analysis/analysis_output/t1d_tpm_vs_apn_w_mei.csv" 
   dbms=csv replace;
run;


