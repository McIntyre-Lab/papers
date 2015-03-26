/* 
 * I am comparing the truncated 2011-07-05/s_3_TTAGGC vs the full dataset to see
 * how coverage counts compare.
 */
libname fru "!MCLAB/arbeitman/arbeitman_fru_network/sasdata";

 data WORK.missing ;
 infile '!MCLAB/arbeitman/arbeitman_fru_network/pipeline_output/coverage_on_fusions/coverage_on_fusions_2011-07-05_3_AH_Male_FruM_A_3_missing.csv'
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=1 ;
    informat fusion_id $9. ;
    informat sample_id $20. ;
    informat mapped_reads best32. ;
    informat reads_in_exon best32. ;
    informat coverage_in_exon best32. ;
    informat exon_length best32. ;
    informat apn best32. ;
    informat rpkm best32. ;
    format fusion_id $9. ;
    format sample_id $20. ;
    format mapped_reads best12. ;
    format reads_in_exon best12. ;
    format coverage_in_exon best12. ;
    format exon_length best12. ;
    format apn best12. ;
    format rpkm best12. ;
 input
             fusion_id $
             sample_id
             mapped_reads
             reads_in_exon
             coverage_in_exon
             exon_length
             apn
             rpkm
 ;
 run;

proc sort data=missing ;
    by apn;
    run;

data fragment;
    set missing;
    rename apn =apn_miss;
    if apn>5 then flag_miss_oops=1;
    else flag_miss_oops=0;
    keep fusion_id apn flag_miss_oops;
    run;

    proc freq data=fragment;
    tables flag_miss_oops;
    run;

    proc univariate data=fragment normal plot;
    var apn_miss;
    run;

proc sort data=fru.results_plus_gov2;
    by fusion_id;
    run;

proc sort data=fragment;
    by fusion_id;
    run;

data check_fragment oops;
    merge fru.results_plus_gov2(in=in1) fragment (in=in2);
    by fusion_id;
    if in1 and not in2 then output oops;
    else output check_fragment;
    run;

proc contents data=check_fragment;
run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_11_20 flag_miss_oops;
    output out=check n=n;
    run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_12_20 flag_miss_oops;
    output out=check n=n;
    run;

********************************************************************************;
data fragment;
    set missing;
    rename apn =apn_miss;
    if apn>10 then flag_miss_oops=1;
    else flag_miss_oops=0;
    keep fusion_id apn flag_miss_oops;
    run;

    proc freq data=fragment;
    tables flag_miss_oops;
    run;

    proc univariate data=fragment normal plot;
    var apn_miss;
    run;

proc sort data=fru.results_plus_gov2;
    by fusion_id;
    run;

proc sort data=fragment;
    by fusion_id;
    run;

data check_fragment oops;
    merge fru.results_plus_gov2(in=in1) fragment (in=in2);
    by fusion_id;
    if in1 and not in2 then output oops;
    else output check_fragment;
    run;

proc contents data=check_fragment;
run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_11_20 flag_miss_oops;
    output out=check n=n;
    run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_12_20 flag_miss_oops;
    output out=check n=n;
    run;

********************************************************************************;
data fragment;
    set missing;
    rename apn =apn_miss;
    if apn>100 then flag_miss_oops=1;
    else flag_miss_oops=0;
    keep fusion_id apn flag_miss_oops;
    run;

    proc freq data=fragment;
    tables flag_miss_oops;
    run;

    proc univariate data=fragment normal plot;
    var apn_miss;
    run;

proc sort data=fru.results_plus_gov2;
    by fusion_id;
    run;

proc sort data=fragment;
    by fusion_id;
    run;

data check_fragment oops;
    merge fru.results_plus_gov2(in=in1) fragment (in=in2);
    by fusion_id;
    if in1 and not in2 then output oops;
    else output check_fragment;
    run;

proc contents data=check_fragment;
run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_11_20 flag_miss_oops;
    output out=check n=n;
    run;

proc means data=check_fragment noprint;
    class flag_fdr_p_contrast_12_20 flag_miss_oops;
    output out=check n=n;
    run;


data panic;
	set check_fragment;
	if apn_miss ge 5;
    keep symbol_cat;
	run;

proc sort data=panic nodupkey;
    by symbol_cat;
    run;

proc print data=panic;
	run;



