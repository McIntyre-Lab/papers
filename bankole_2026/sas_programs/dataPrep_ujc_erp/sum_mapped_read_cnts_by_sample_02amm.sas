

/*

mapped reads per sample
reads per sample thru ujc
  
*/



libname sex "!MCLAB/sex_specific_splicing/sasdata";
libname local "/TB20/sex_specific_splicing/sasdata";



/* import number mapped reads */
%macro imp_mapped (dir, input, species) ;

proc import datafile = "!MCLAB/&dir./qc_output/&input."
out = &species._mapped
dbms = csv replace ;
guessingrows = MAX ;
run;

data &species._mapped2 ;
set &species._mapped ;
species = scan(sample, 1, '_') ;
sex = scan(sample, 2, '_') ;
rep = scan(sample, 3, '_') ;
techRep = scan(sample, 4, '_') ;
if find(sample, "dser") ge 1 then newID = compress(species||'_'||sex||'_'||rep) ;
else newID = compress("d"||species||'_'||sex||'_'||rep) ;
drop sample ;
rename newID = sample ;
run ;

/* sum to sample gene */
proc sort data = &species._mapped2;
by sample ;
run;

proc means data = &species._mapped2;
by sample ;
var mapped_read_num input_read_num ;
output out = &species._sum_mapped sum =/ autoname;
run;

data &species._mapped_read_cnts ;
set &species._sum_mapped ;
drop _type_ _freq_ ;
run;

proc export data = &species._mapped_read_cnts 
outfile = "!MCLAB/&dir./qc_output/&species._mapped_read_cnts_per_sample.csv"
dbms = csv replace ;
run ;
%mend ;

%imp_mapped (lmm_rmg_dros_data, mel_data_2_dmel6_mapped_read_cnts.csv, dmel) ;
%imp_mapped (lmm_rmg_dros_data, sim_data_2_dsim2_mapped_read_cnts.csv, dsim) ;
%imp_mapped (lmm_axk_head_data, dser_data_2_dser1_mapped_read_cnts.csv, dser) ;


/* import number mapped reads  for yak and san! */
%macro imp_mapped (dir, input, species) ;

proc import datafile = "!MCLAB/&dir./qc_output/&input."
out = &species._mapped
dbms = csv replace ;
guessingrows = MAX ;
run;

data &species._mapped2 ;
set &species._mapped ;
species = scan(sample, 1, '_') ;
sex = scan(sample, 2, '_') ;
rep = "rep1" ;
techRep = scan(sample, 3, '_') ;
newID = compress(species||'_'||sex||'_'||rep) ;
drop sample ;
rename newID = sample ;
run ;

/* sum to sample gene */
proc sort data = &species._mapped2;
by sample ;
run;

proc means data = &species._mapped2;
by sample ;
var mapped_read_num input_read_num ;
output out = &species._sum_mapped sum =/ autoname;
run;

data &species._mapped_read_cnts ;
set &species._sum_mapped ;
drop _type_ _freq_ ;
run;

proc export data = &species._mapped_read_cnts 
outfile = "!MCLAB/&dir./qc_output/&species._mapped_read_cnts_per_sample.csv"
dbms = csv replace ;
run ;
%mend ;

%imp_mapped (lmm_rlr_head_data, dyak_data_2_dyak2_mapped_read_cnts.csv, dyak) ;
%imp_mapped (lmm_rlr_head_data, dsan_data_2_dsan1_mapped_read_cnts.csv, dsan) ;


/* set all species together */
data all_mapped ;
set dmel_mapped_read_cnts dsim_mapped_read_cnts dser_mapped_read_cnts dsan_mapped_read_cnts dyak_mapped_read_cnts ;
run;

proc freq data = all_mapped ;
tables sample / out = sample_cnt ;
run ;

/* count reads per sample from ujc count file - AFTER id_ujc */
%macro cnt_check (species, genome) ;

proc sort data = sex.&species._ujc_cnts_sumtr_stack ;
by sample ;
run ;

proc means data = sex.&species._ujc_cnts_sumtr_stack;
by sample ;
var jxnHash_sumReadNum;
output out = &species._ujc_sum_sample sum =/ autoname;
run;

data &species._ujc_sum_sample2 ;
set &species._ujc_sum_sample ;
drop _type_ _freq_ ;
run;
%mend ;

%cnt_check (dmel, dmel6) ;
%cnt_check (dsim, dsim2) ;
%cnt_check (dser, dser1) ;
%cnt_check (dsan, dsan1) ;
%cnt_check (dyak, dyak2) ;


data all_ujc_sums ;
set dmel_ujc_sum_sample2 dsim_ujc_sum_sample2 dser_ujc_sum_sample2 dsan_ujc_sum_sample2 dyak_ujc_sum_sample2 ;
run;

proc sort data = all_mapped ;
by sample ;
proc sort data = all_ujc_sums ;
by sample ;
run;

data mapped_jxnHash_readCnts ;
merge all_mapped (in=in1) all_ujc_sums (in=in2) ;
by sample ;
run; 

data sex.mapped_jxnHash_readCnt_summary ;
retain sample input_read_num_sum mapped_read_num_Sum prop_mapped jxnHash_sumReadNum_sum  prop_in_jxnHash;
set mapped_jxnHash_readCnts ;
prop_mapped = mapped_read_num_Sum/input_read_num_sum ;
prop_in_jxnHash = jxnHash_sumReadNum_sum/mapped_read_num_Sum ;
run ;

proc export data = sex.mapped_jxnHash_readCnt_summary 
outfile = "!MCLAB/sex_specific_splicing/mapped_jxnHash_readCnt_summary.csv"
dbms = csv replace ;
run;

