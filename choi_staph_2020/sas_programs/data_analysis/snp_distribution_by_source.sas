
libname snp 'Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata';
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";


	data match;
	set snp.relapse_snp_counts_match;

	data no_match;
	set snp.relapse_snp_counts_nomatch;
	run;



	data pairs;
	set match no_match;
	isolate1=substr(pair,6,4);
	isolate2=substr(pair,11,4);
	id1=isolate1*1;
	id2=isolate2*1;
	drop cc st;
	run;


data choi_841;
set choi.choi_841;
keep id ridom_cc female age episode_number st patientnumber;
run;

data id1;
set choi_841;
rename id=id1 patientnumber=pn1 ridom_cc=cc1 female=female1 age= age1 episode_number=epi1 st=st1;
run;

data id2;
set choi_841;
rename id=id2 patientnumber=pn2 ridom_cc=cc2 female=female2 age= age2 episode_number=epi2 st=st2;
run;

proc sort data=pairs ;
by id1;

proc sort data=id1;
by id1;

data pairs_id1;
merge pairs(in=in1) id1(in=in2);
by id1;
if in1;

proc sort data=pairs_id1;
by id2;

proc sort data=id2;
by id2;

data pairs_id1_id2;
merge pairs_id1 (in=in1) id2(in=in2);
by id2;
if in1;
if sum_no_match le 100 then flag_le_100=1;
else flag_le_100=0;
if st1=st2 then flag_st_match=1;
else flag_st_match=0;
run;


proc freq data =pairs_id1_id2;
tables flag_st_match*flag_le_100;
tables flag_le_100*flag_is_a_pair;
run;

data check;
set pairs_id1_id2;
where flag_is_a_pair='0' and flag_st_match=1 and flag_le_100=1;
run;


proc freq data=no_match1;
tables ref;
run;

proc means data=pairs_id1_id2 min p10 median p90 max;
class flag_st_match flag_is_a_pair;
var sum_no_match;
run;

/*note that I checked that all with le 100 snps have same st type*/

data for_plotting_snps;
set pairs_id1_id2;
if flag_is_a_pair='1' then group="Pair";
else if flag_le_100=1 then group="PFGE";
else if flag_st_match=1 then group="ST";
else group="CC";
keep sum_no_match group ref;
run;

data snp_pair;
set for_plotting_snps;
where flag_is_a_pair='1';
run;

data snp_st_match;
set for_plotting_snps;
where 
run;

data snp_cc_match;
set for_plotting_snps;
where flag_st_match=0 ;
run;


PROC EXPORT DATA= WORK.for_plotting_snps
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\for_plotting_snps.csv" 
            DBMS=CSV REPLACE;
RUN;



PROC EXPORT DATA= WORK.match1
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\match.csv" 
            DBMS=CSV REPLACE;
RUN;
