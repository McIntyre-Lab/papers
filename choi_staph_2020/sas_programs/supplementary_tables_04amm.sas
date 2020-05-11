libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
libname TB14 "/home/ammorse/TB14/staph_relapse/sasdata";
libname lmm "!MCLAB/staph/seoung_ho_project/sasdata/data_analysis";

/*
libname relapse "/home/ammorse/staph_relapse/sasdata";
libname lmm  "/home/ammorse/staph_relapse/sasdata/data_analysis";
*/
/*

inputs:
    lmm.id_link_2_patient_pfge
    relapse.sampleID_st_cc_ref

outputs:
    relapse.sampleID_st_cc_ref_kitchen_sink
    relapse.sampleID_st_cc_ref_4_sra

    !MCLAB/staph/seoung_ho_project/supplementary_material/sampleID_st_cc_ref_4_sra.csv
*/



/* table containing new sra_sampleID, ST, CC, reference info */

data lmm_table ;
set lmm.id_link_2_patient_pfge ;
rename id = sampleID ;
drop num_pairs ;
run;

proc freq data = lmm_table ;
table patientnumber / out = patientcnts ;
run;  /* 69 patient numbers starting with 0.... */

data prep_table ;
set relapse.sampleID_st_cc_ref ;
*where pairNum ne 5 ;
drop  flag_ref ;
run ;

proc sort data = prep_table ;
by pairNum ;
run ;

/* merge lmm table with my table */
proc sort data = lmm_table ;
by sampleID ;
proc sort data = prep_table ;
by sampleID ;
run;

data relapse_table ;
merge lmm_table (in=in1) prep_table (in=in2) ;
by sampleID ;
if in2 ;
run;  

proc sort data = relapse_table ;
by patientnumber ;
run;

/* 0 in oops, 76 patientnumber with no seq data */
proc freq data = relapse_table ;
table reference sampleID patientnumber;
run;
/*
reference                               Frequency     
-------------------------------------------------------
CA_347                                         2
ED98                                          24
MSSA476                                        4
Newman                                        40
ST20130941                                     6
TCH60                                          2

                                
sampleIDs unique   
patient number not unique (expected) */

proc sort data = relapse_table  ;
by  pairNum ;
run;

data relapse_table2 ;
set relapse_table ;
	count +1 ;
	by pairNum ;
	if first.pairNum then count = 1;
run ;

data relapse.sampleID_st_cc_ref_kitchen_sink;
retain sra_sampleID patientnumber sampleID sampleNum ST CC reference taxonomy_ID  ;
set relapse_table2 ;
if count = 1 then sra_sampleID = compress(patientnumber||'A') ;
else if count = 2 then sra_sampleID = compress(patientnumber||'B') ;
rename taxonomy_ID = ref_taxonomyID  ;
if pairNum = 5 and flag_failed_library = 0 then flag_failed_library_pair = 1 ;
    else flag_failed_library_pair = 0;
drop count ;
run ;

data relapse.sampleID_st_cc_ref_4_sra ;
set  relapse.sampleID_st_cc_ref_kitchen_sink;
drop sampleID sampleNum pairNum pfge;
where flag_failed_library = 0 ;
run ;

proc freq data = relapse.sampleID_st_cc_ref_4_sra  ;
tables ST flag_failed_library ;
run;  /*
   ST    Frequency
 -----------------
    1           4
    5          14
    8          36
   15           6
   30           2
   45           2
  105           9
  840           2
 1181           4    */

proc freq data = relapse.sampleID_st_cc_ref_4_sra  ;
tables CC ;
run;  /*
CC    Frequency
---------------
 1           4
 5          25
 8          40
15           6
30           2
45           2  */

proc sort data = relapse.sampleID_st_cc_ref_4_sra ;
by  sra_sampleID ;
run;


proc export data = relapse.sampleID_st_cc_ref_4_sra 
/*outfile = "/home/ammorse/staph_relapse/supplementary_material/sampleID_st_cc_ref_4_sra.csv"*/
outfile = "!MCLAB/staph/seoung_ho_project/supplementary_material/sampleID_st_cc_ref_4_sra.csv" 
dbms = csv replace ;
run ;




