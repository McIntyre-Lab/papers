
libname snp 'Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata';
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";


	data match;
	set snp.relapse_snp_counts_match;
	drop cc st pairnum;

	data no_match;
	set snp.relapse_snp_counts_nomatch;
	drop st cc pairnum;
	run;


data pair;
set match no_match;
id1 = input(pair1, 12.);
id2= input(pair2,12.);
drop pair1 pair2;
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

proc sort data=pair;
by id1;

proc sort data=id1;
by id1;

data pairs_id1;
merge pairs(in=in1) id1(in=in2);
by id1;
if in1 and in2;
run;


proc sort data=pairs_id1;
by id2;

proc sort data=id2;
by id2;

data pairs_id1_id2;
merge pairs_id1 (in=in1) id2(in=in2);
by id2;
if in1 and in2;
if sum_no_match le 100 then flag_le_100=1;
else flag_le_100=0;
if st1=st2 then flag_st_match=1;
else flag_st_match=0;
run;

proc means data=pairs_id1_id2 min p10 median p90 max;
class flag_st_match flag_is_a_pair;
var sum_no_match;
run;

/*note that I checked that all with le 100 snps have same st type*/

data snp.for_plotting_snps;
set pairs_id1_id2;
if flag_is_a_pair='1' then group="Pair";
else if flag_le_100=1 then group="PFGE";
else if flag_st_match=1 then group="ST";
else group="CC";
keep sum_no_match group ref;
run;


PROC EXPORT DATA= WORK.for_plotting_snps
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\for_plotting_snps1.csv" 
            DBMS=CSV REPLACE;
RUN;


