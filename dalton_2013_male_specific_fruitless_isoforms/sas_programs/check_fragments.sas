libname fru "Z:\arbeitman\arbeitman_fru_network\sasdata";


data fragment;
set fru.missing;
rename apn =apn_miss;
if apn>5 then flag_miss_oops=1;
else flag_miss_oops=0;
keep fusion_id apn flag_miss_oops;
run;

proc freq data=fragment;
tables flag_miss_oops;
run;


proc univariate data=fragment normal plot;
var apn;
run;


proc sort data=fru.results_plus_gov2;
by fusion_id;

proc sort data=fragment;
by fusion_id;

data check_fragment oops;
merge fru.results_plus_gov2(in=in1) fragment (in=in2);
by fusion_id;
if in1 and not in2 then output oops;
else output check_fragment;
run;

ods results;
proc contents data=check_fragment;
run;

proc means data=check_fragment noprint;
class
    flag_fdr_p_contrast_10_20
	flag fdr_p_contrast_11_20
	flag_miss_oops;
	output out=check n=n;
	run;

	data panic;
	set check_fragment;
	if apn_miss ge 5;
	run;


	proc print data=panic;
	var gene;
	run;


