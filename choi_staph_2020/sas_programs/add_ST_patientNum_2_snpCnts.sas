

/*
add patient numbe,  ST etc to snp counts 
    match and noMatch


*/



/* add patient number for each pair to matching pairs */
data add_pn ;
set relapse.sampleID_st_cc_ref_kitchen_sink;
pair1 = put(sampleID, 4.) ;
drop sampleID ;
where flag_failed_library_pair = 0 and flag_failed_library =  0;
keep pairNum pair1  patientnumber CC ST ref_taxonomyID reference ;
run ;

proc sort data = add_pn nodups;
by _all_ ;
proc sort data = add_pn;
by pair1 ;
run;

data iso1 ;
retain pair1 pair2 ;
set relapse.relapse_snp_counts_match ;
run ;

proc sort data = iso1  ;
by pair1 ;
run ;

data add ;
merge iso1  (in=in1) add_pn (in=in2) ;
by pair1 ;
if in1 ;
run ;

data add2 ;
set add ;
rename ST = pair1_ST ;
run;

data add_pn_2 ;
set add_pn ;
rename pair1 = pair2 ;
run ;


proc sort data = add2 ;
by pair2 ;
proc sort data = add_pn_2 ;
by pair2 ;
run ;

data add3 ;
merge add2 (in=in1) add_pn_2 (in=in2) ;
by pair2 ;
if in1 ;
run;

data add4 ;
retain pair1 pair2 pair pair1_ST ST ;
set add3 ;
rename ST = pair2_ST ;
run;


data ready_match ;
retain pair1 pair2 pair patientnumber pair1_ST pair2_ST CC reference ref_taxonomyID patientnumber sra_pair ;
set add4 ;
sra_pair = compress("pair_"||patientnumber||"A_"||patientnumber||"B") ;
rename sum_match = sum_snp_match ;
rename sum_no_match = sum_snps_no_match ;
rename total = total_snp ;
drop ref pairNum  ;
run;



/* add ST to no Match */
data add_pn ;
set relapse.sampleID_st_cc_ref_kitchen_sink;
pair1 = put(sampleID, 4.) ;
drop sampleID ;
where flag_failed_library_pair = 0 and flag_failed_library =  0;
keep pair1  CC ST  ;
run ;

proc sort data = add_pn nodups;
by _all_ ;
proc sort data = add_pn;
by pair1 ;
run;

data iso_noM ;
retain pair1 pair2 ;
set relapse.relapse_snp_counts_noMatch ;
run ;

proc sort data = iso_noM ;
by pair1 ;
run ;

data add_noM ;
merge iso_noM  (in=in1) add_pn (in=in2) ;
by pair1 ;
if in1 ;
run ;

data add2_noM ;
set add_noM ;
rename ST = pair1_ST ;
run;

data add_pn_2 ;
set add_pn ;
rename pair1 = pair2 ;
run ;


proc sort data = add2_noM ;
by pair2 ;
proc sort data = add_pn_2 ;
by pair2 ;
run ;

data add3_noM ;
merge add2_noM (in=in1) add_pn_2 (in=in2) ;
by pair2 ;
if in1 ;
run;

data add4_noM ;
retain pair1 pair2 pair pair1_ST ST ;
set add3_noM ;
rename ST = pair2_ST ;
run;


data ready_noMatch ;
retain pair1 pair2 pair pair1_ST pair2_ST CC ref ;
set add4 ;
rename ref = reference ;
rename sum_no_match = sum_snps_no_match ;
rename total = total_snp ;
run;


/* make perm */
data relapse.snp_cnts_match_ST_CC_ref_sra ;
retain pair1 pair2 pair sra_pair patientNumber pair1_st pair2_st CC reference ;
set ready_match ;
run;

data relapse.snp_cnts_noMatch_ST_CC_ref ;
set ready_noMatch ;
run;


proc export data = relapse.snp_cnts_match_ST_CC_ref_sra 
/*outfile = "/home/ammorse/staph_relapse/supplementary_material/relapse_snp_cnts_ST_CC_ref_sra.csv"*/
outfile = "!MCLAB/staph/seoung_ho_project/supplementary_material/relapse_snp_cnts_match_ST_CC_ref_sra.tsv" 
dbms = csv replace ;
run ;


