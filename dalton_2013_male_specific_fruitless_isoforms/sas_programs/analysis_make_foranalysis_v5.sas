/*
 * REVISIONS: 12/21/2011
 *            - Changed libname
 *            - Added Gene information from demel 5.30 for use in analysis model
 */

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

data all_coverage_w_flag;
    set fru.all_coverage_counts_with_key;
    if index(trt, 'Fru') > 0 then flag_fru = 1; else flag_fru = 0;
    if index(trt, 'dsx') > 0 then flag_dsx = 1; else flag_dsx = 0;
    if index(trt, 'CS') > 0 or index(trt,'Ber') > 0 then flag_control = 1; else flag_control = 0;
    run;
        
proc sort data=all_coverage_w_flag;
    by fusion_id sample_id;
    run;

proc sort data=fru.flag_no_trt;
    by fusion_id sample_id;
    run;

data short_list oops;
    merge all_coverage_w_flag (in=in1) fru.flag_no_trt (in=in2);
    by fusion_id sample_id;
    if in1 and in2 then output short_list;
    else output oops;
    run;

data fru;
    set short_list;
    if flag_no_trt = 0 and flag_dsx = 0;
    drop flag_no_trt flag_dsx;
    run;

data dsx;
    set short_list;
    where flag_no_trt = 0 and flag_fru = 0;
    drop flag_no_trt flag_fru;
    run;

proc sort data=fru;
    by fusion_id;
    run;

proc sort data=dsx;
    by fusion_id;
    run;

proc sort data=fru.flag_no_var;
    by fusion_id;
    run;

proc sort data=fru.flag_zero_mean;
    by fusion_id;
    run;

data fusion2gene;
    set dmel530.Fb530_si_fusions_unique;
    if Genes_per_fusion > 1 then multigene = 1;
    else multigene = 0;
    keep fusion_id symbol_cat multigene;
    run;

proc sort data=fusion2gene;
    by fusion_id;
    run;

data foranalysis_fru;
    merge fru (in=in1) 
          fru.flag_no_var (in=in2) 
          fru.flag_zero_mean (in=in3)
          fusion2gene (in=in4);
    by fusion_id;
    if in1 and in2 and in3 and in4;
    run; * 2291058 obs;

data foranalysis_dsx;
    merge dsx (in=in1) 
          fru.flag_no_var (in=in2) 
          fru.flag_zero_mean (in=in3)
          fusion2gene (in=in4);
    by fusion_id;
    if in1 and in2 and in3 and in4;
    run; * 2291058 obs;

data fru.foranalysis_no_multi_fru;
    retain symbol fusion_id;
    set foranalysis_fru;
    where multigene = 0 and flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1722882 obs;

data fru.foranalysis_with_multi_fru;
    retain symbol fusion_id;
    set foranalysis_fru;
    where  flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1801238 obs;

data fru.foranalysis_no_multi_dsx;
    retain symbol fusion_id;
    set foranalysis_dsx;
    where multigene = 0 and flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1722882 obs;

data fru.foranalysis_with_multi_dsx;
    retain symbol fusion_id;
    set foranalysis_dsx;
    where  flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1801238 obs;
