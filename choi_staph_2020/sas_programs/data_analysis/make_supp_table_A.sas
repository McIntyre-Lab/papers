
libname snp 'Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata';
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";

libname snp "!MCLAB/staph/seoung_ho_project/sasdata";



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
keep id ridom_cc st patientnumber;
run;

data id1;
set choi_841;
rename id=id1 patientnumber=pn1 ridom_cc=cc1 st=st1;
run;

data id2;
set choi_841;
rename id=id2 patientnumber=pn2 ridom_cc=cc2 st=st2;
run;

proc sort data=pair;
by id1;

proc sort data=id1;
by id1;

data pairs_id1;
merge pair(in=in1) id1(in=in2);
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
rename ref=reference;
run;

proc sort data=snp.cc_refs;
by reference;

proc sort data=pairs_id1_id2;
by reference;

data pairs_with_ref_info;
merge pairs_id1_id2 snp.cc_refs;
by reference;

rename CC=reference_CC ST=Reference_ST genotype_pair_0_0=Both_Match_Reference
genotype_pair_1_1=Both_SNP cc1=SpaType1 cc2=SpaType2 
st1=MLST_ST1 ST2=MLST_ST2;

drop genotype_pair_0_1 genotype_pair_1_0 
 pair match_freq nomatch_freq total
id1 id2 flag_ref;

run;

data snp.supp_table_A;
retain reference taxonomy_ID reference_CC reference_ST flag_is_a_pair 
pn1 pn2   MLST_ST1 MLST_ST2 SpaType1 Spatype2 
Both_Match_Reference Both_SNP  sum_match sum_no_match;
set pairs_with_ref_info;
run;

proc sort data=snp.supp_table_A;
by sum_no_match;
run;

proc export data = snp.supp_table_A
outfile = "!MCLAB/staph/seoung_ho_project/supplementary_material/supplementary_table_A.csv"
dbms = csv replace ;
run;


