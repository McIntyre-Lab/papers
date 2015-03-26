
libname fru '!MCLAB/Fru_network/sasdata';

/* create dataset to make means, remove samples without trt information */
data formeans;
    set fru.all_coverage_counts_with_key;
    where trt ne "";
    run;

proc sort data=formeans;
    by fusion_id trt;
    run;

/* Calculate Mean and SD for LogRPKM+100 and transpose */
proc means data=formeans noprint;
    var logrpkm;
    by fusion_id;
    class trt;
    output out=logrpkm_means mean=mean stddev=SD;
    run;

data logrpkm_means;
    set logrpkm_means;
    where _type_ = 1;
    run;

proc transpose data=logrpkm_means out=logtrans_means suffix=_mean_logrpkm;
    by fusion_id;
    var mean;
    id trt;
    run;

proc transpose data=logrpkm_means out=logtrans_sd suffix=_SD_logrpkm;
    by fusion_id;
    var SD;
    id trt;
    run;

/* Calculate Mean and SD for RPKM and transpose */
proc means data=formeans noprint;
    var rpkm;
    by fusion_id;
    class trt;
    output out=rpkm_means mean=mean stddev=SD;
    run;

data rpkm_means;
    set rpkm_means;
    where _type_ = 1;
    run;

proc transpose data=rpkm_means out=rpkmtrans_means suffix=_mean_rpkm;
    by fusion_id;
    var mean;
    id trt;
    run;

proc transpose data=rpkm_means out=rpkmtrans_sd suffix=_SD_rpkm;
    by fusion_id;
    var SD;
    id trt;
    run;

/* Combine logrpkm+100 and rpkm mean and SD into one dataset */
proc sort data=logtrans_means;
    by fusion_id;
    run;

proc sort data=logtrans_sd;
    by fusion_id;
    run;

proc sort data=rpkmtrans_means;
    by fusion_id;
    run;

proc sort data=rpkmtrans_sd;
    by fusion_id;
    run;

data fru.all_meansv2;
    retain fusion_id
           AH_BerF_mean_logrpkm
           AH_BerF_SD_logrpkm
           AH_BerM_mean_logrpkm
           AH_BerM_SD_logrpkm
           AH_CS_mean_logrpkm
           AH_CS_SD_logrpkm
           AH_CSFemale_mean_logrpkm
           AH_CSFemale_SD_logrpkm
           AH_Female_FruM_A__mean_logrpkm
           AH_Female_FruM_A__SD_logrpkm
           AH_Female_FruM_B__mean_logrpkm
           AH_Female_FruM_B__SD_logrpkm
           AH_Female_FruM_C__mean_logrpkm
           AH_Female_FruM_C__SD_logrpkm
           AH_FruP14_440_mean_logrpkm
           AH_FruP14_440_SD_logrpkm
           AH_FruW12_ChaM5_mean_logrpkm
           AH_FruW12_ChaM5_SD_logrpkm
           AH_Male_FruM_A__mean_logrpkm
           AH_Male_FruM_A__SD_logrpkm
           AH_Male_FruM_B__mean_logrpkm
           AH_Male_FruM_B__SD_logrpkm
           AH_Male_FruM_C__mean_logrpkm
           AH_Male_FruM_C__SD_logrpkm
           AH_dsxD_mean_logrpkm
           AH_dsxD_SD_logrpkm
           AH_dsxNullF_mean_logrpkm
           AH_dsxNullF_SD_logrpkm
           AH_dsxNullM_mean_logrpkm
           AH_dsxNullM_SD_logrpkm
           AH_BerF_mean_rpkm
           AH_BerF_SD_rpkm
           AH_BerM_mean_rpkm
           AH_BerM_SD_rpkm
           AH_CS_mean_rpkm
           AH_CS_SD_rpkm
           AH_CSFemale_mean_rpkm
           AH_CSFemale_SD_rpkm
           AH_Female_FruM_A__mean_rpkm
           AH_Female_FruM_A__SD_rpkm
           AH_Female_FruM_B__mean_rpkm
           AH_Female_FruM_B__SD_rpkm
           AH_Female_FruM_C__mean_rpkm
           AH_Female_FruM_C__SD_rpkm
           AH_FruP14_440_mean_rpkm
           AH_FruP14_440_SD_rpkm
           AH_FruW12_ChaM5_mean_rpkm
           AH_FruW12_ChaM5_SD_rpkm
           AH_Male_FruM_A__mean_rpkm
           AH_Male_FruM_A__SD_rpkm
           AH_Male_FruM_B__mean_rpkm
           AH_Male_FruM_B__SD_rpkm
           AH_Male_FruM_C__mean_rpkm
           AH_Male_FruM_C__SD_rpkm
           AH_dsxD_mean_rpkm
           AH_dsxD_SD_rpkm
           AH_dsxNullF_mean_rpkm
           AH_dsxNullF_SD_rpkm
           AH_dsxNullM_mean_rpkm
           AH_dsxNullM_SD_rpkm;
    merge logtrans_means
          logtrans_sd
          rpkmtrans_means
          rpkmtrans_sd;
    by fusion_id;
    drop _name_;
    run;


