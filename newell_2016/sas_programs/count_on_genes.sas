

libname ribo 'z:\arbeitman\arbeitman_ribotag\sas_data';
libname dmel 'z:\useful_dmel_data\flybase551\sasdata';

libname ribo '/home/mcintyre/S/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/sas_data';

proc contents data=ribo.On_calls_gt_apn0;
run;

proc freq data=ribo.On_calls_gt_apn0;
 tables flag_IPfemale_on   
 flag_IPmale_on    
 flag_InputFemale_on    
 flag_InputMale_on     
 flag_fusion_all_on0    
 flag_fusion_on0 
;
run;

proc sort data=ribo.On_calls_gt_apn0;
by fusion_id;

proc contents data=dmel.fb551_si_fusions_unique_flagged;
  run;

proc sort data=dmel.fb551_si_fusions_unique_flagged;
  by fusion_id;
  run;

  data count_on_genes oops;
  merge ribo.On_calls_gt_apn0 (in=in1) dmel.fb551_si_fusions_unique_flagged (in=in2);
by fusion_id;
if in1 and in2 then output count_on_genes;
else output oops;
run;

*there are 63706 fusions in the dmel.fb551_si_fusions_unique_flagged  included in this file are redundant fusions that have been eliminated. There are only 63181 non-redundant fusions;

proc contents data=count_on_genes;
run;

data no_multi;
set count_on_genes;
where Genes_per_fusion=1;
run;
*6030 fusions removed;

data all_exons;
set no_multi;
keep symbol_cat;

proc sort data=all_exons out=all_genes nodupkey;
by symbol_cat;
run;


data all_exons_analyzed;
set ribo.results_by_fusion_new_model;
*from gene count import;
keep fusion_id symbol_cat;

proc sort data=all_exons_analyzed out=all_genes_analyzed nodupkey;
by symbol_cat;
run;


*14716 genes in the database ;

%macro on_gene(trt);

	proc freq data=no_multi noprint;
	tables flag_&trt._on  *symbol_cat/out=count_&trt._on;
	run;

	data &trt._on_genes;
	set count_&trt._on;
	where flag_&trt._on=1;
	keep symbol_cat flag_&trt._on;
	run;

	proc sort data=&trt._on_genes;
	by symbol_cat;
	run;
%mend;

%on_gene(inputfemale);
%on_gene(inputmale);
%on_gene(IPfemale);
%on_gene(IPmale);


data count_on_genes2;
merge IPfemale_on_genes(in=in1) IPmale_on_genes(in=in2)
inputfemale_on_genes(in=in3) inputmale_on_genes(in=in4) all_genes (in=in5)
all_genes_analyzed (in=in6);
by symbol_cat;
if in1 and in2 and in3 and in4 then gene_on="all_on";
	else if in1 or in2 or in3 or in4 then gene_on="some";
	else gene_on="none";
if in1 or in2 then flag_trap_on=1;
	else flag_trap_on=0;
if in3 or in4 then flag_input_on=1;
	else flag_input_on=0;
if flag_inputfemale_on="." then flag_inputfemale_on=0;
if flag_inputmale_on="." then flag_inputmale_on=0;
if flag_ipfemale_on="." then flag_ipfemale_on=0;
if flag_ipmale_on="." then flag_ipmale_on=0;

if flag_inputfemale_on=1 or flag_inputmale_on=1 or
flag_IPfemale_on=1 or flag_IPmale_on=1 then flag_gene_detected=1;
	else flag_gene_detected=0;
	if in6 then flag_analyzed=1;
	else flag_analyzed=0;
run;

proc freq data=count_on_genes2;
tables gene_on*flag_gene_detected;
tables flag_gene_detected*flag_trap_on;
tables flag_gene_detected*flag_input_on;
run;

proc freq data=count_on_genes2;
*where flag_gene_detected=1 and flag_input_on=0;
tables flag_IPfemale_on*flag_IPmale_on;
run;



proc freq data=count_on_genes2;
*where flag_gene_detected=1 and flag_input_on=0;
tables flag_Inputfemale_on*flag_Inputmale_on;
run;
proc freq data=count_on_genes2;
where gene_on="some";
tables flag_trap_on*flag_input_on;
run;

proc print data=count_on_genes2;
where flag_trap_on=1 and flag_input_on=0;
var symbol_cat flag_ipfemale_on flag_ipmale_on;
run;

proc freq data= count_on_genes2;
*where flag_input_on=1;
where gene_on="all_on";
tables flag_ipfemale_on*flag_ipmale_on;
run;

data check1;
set ribo.results_by_fusion_new_model;
where symbol_cat ="fru";
run;


data check2;
set all_genes_analyzed;
where symbol_cat ="fru";
run;

proc freq data= count_on_genes2;
tables gene_on*flag_analyzed;
run;


data check;
set  count_on_genes2;
where flag_analyzed=0 and gene_on="all_on";
run;

*look at exons that are detected in maleip but not female ip;
data detect_maleip;
set count_on_genes2;
where flag_Ipfemale_on=0 and flag_Ipmale_on=1 ;
run;

proc freq data= detect_maleip;
tables flag_Inputfemale_on*flag_Inputmale_on;
run;

proc freq data= detect_maleip noprint;
tables symbol_cat/out=checking_male_detect;
run;
*3293  genes!;

proc sort data=all_results;
by symbol_cat;
proc sort data=detect_maleip;
by symbol_cat;

data check_male_list1 oops1 oops2;
merge detect_maleip(in=in1) all_results(in=in2);
by symbol_cat;
if in1 and in2 then output oops1;
else if in1 then output check_male_list1;
else output oops2;
run;



data pattern;
set ribo.on_calls_gt_apn0;
pattern=cat(trim(flag_Inputfemale_on)||trim(flag_inputmale_on)||trim(flag_ipfemale_on)||trim(flag_ipmale_on));
run;

proc print data=pattern (obs=10);
run;

proc freq data=pattern;
tables pattern;

proc sort data=pattern;
by fusion_id;

proc sort data=all_exons_analyzed;
by fusion_id;

data count_genes_by_pattern;
merge pattern (in=in1) all_exons_analyzed;
by fusion_id;
if in1;
run;

proc freq data=count_genes_by_pattern noprint;
tables symbol_cat*pattern/out=count_both;
run;

data count_both1;
set count_both;
rename count=num_exons;
drop percent;
run;

proc freq data=count_both noprint;
tables symbol_cat/out=gene_categories;


proc print data=count_both1(obs=10);
run;

proc sort data=gene_categories;
by symbol_cat;

proc sort data=count_both1;
by symbol_cat;

data how_many_patterns;
merge count_both1 gene_categories;
by symbol_cat;


proc print data=how_many_patterns(obs=10);
run;

proc freq data=how_many_patterns;
where count=1;
tables pattern;
run;

data gene_list;
set how_many_patterns;
where count=1;
drop count percent;
run;

proc sort data=gene_list;
by pattern;
run;

proc export data=gene_list
dbms=csv
outfile='/home/mcintyre/S/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/gene_pattern_unique.csv';
run;





