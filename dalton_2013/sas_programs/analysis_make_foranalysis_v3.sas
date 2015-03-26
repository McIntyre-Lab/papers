/*
 * REVISIONS: 12/21/2011
 *            - Changed libname
 *            - Added Gene information from demel 5.30 for use in analysis model
 */

libname fru '!MCLAB/Fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

proc sort data=short_list2;
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

data foranalysis;
    merge short_list2 (in=in1) 
          fru.flag_no_var (in=in2) 
          fru.flag_zero_mean (in=in3)
          fusion2gene (in=in4);
    by fusion_id;
    if in1 and in2 and in3 and in4;
    run; * 2291058 obs;

data fru.foranalysis_no_multi;
    retain symbol fusion_id;
    set foranalysis;
    where multigene = 0 and flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1722882 obs;

data fru.foranalysis_with_multi;
    retain symbol fusion_id;
    set foranalysis;
    where  flag_no_var =0 and flag_zero_mean = 0;
    drop flag_no_var flag_zero_mean multigene;
    run; *1801238 obs;
