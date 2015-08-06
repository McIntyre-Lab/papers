/*******************************************************************************
* Filename: misc_test_splice_test.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Run a basic splice model to determine if there is a splicing
* effect between treatments
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Compare CS control vs dsxNull */
data dsx;
    set SEM.dsxnull_inr_fusions;
    if trt eq 'AH_BerF' then delete;
    run;

goptions reset=all;
ods listing close;
ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/splicing/dsxNull_splicing_model/InR_cs_vs_null.pdf';
proc glimmix data=dsx plots=studentpanel;
    class fusion_id trt ;
    model logrpkm = trt|fusion_id /htype=1;
    lsmeans trt*fusion_id /slice=trt slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;
ods pdf close;
ods listing;

/* Compare berlin control vs dsxNull */
data dsx;
    set SEM.dsxnull_inr_fusions;
    if trt eq 'AH_CSFemale' then delete;
    run;

goptions reset=all;
ods listing close;
ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/splicing/dsxNull_splicing_model/InR_ber_vs_null.pdf';
proc glimmix data=dsx plots=studentpanel;
    class fusion_id trt ;
    model logrpkm = trt|fusion_id /htype=1;
    lsmeans trt*fusion_id /slice=trt slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;
ods pdf close;
ods listing;

/* Compare combined control vs dsxNull */
data dsx;
    set SEM.dsxnull_inr_fusions;
    run;

goptions reset=all;
ods listing close;
ods pdf file='!MCLAB/cegs_sem_sd_paper/analysis_output/splicing/dsxNull_splicing_model/InR_ctrl_vs_null.pdf';
proc glimmix data=dsx plots=studentpanel;
    class fusion_id trt2 ;
    model logrpkm = trt2|fusion_id /htype=1;
    lsmeans trt2*fusion_id /slice=trt2 slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;
ods pdf close;
ods listing;

