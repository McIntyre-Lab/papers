

libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
filename mymacros "!MCLAB/maize_ozone_final/2014/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);

libname mscramm "!MCLAB/staph/staph_DNA_v2/mscramms/sasdata";
libname ref "!MCLAB/staph/staph_DNA_v2/Nterm_manuscript/analysis_references/sasdata";


/* 
merge in CC information based on ST 

    use to choose reference for snp calling
*/

/* get ST to CC from reference table */

data CC ;
set ref.reference_table_locked ;
keep ST CC ;
if ST = . then delete ;
if CC = . then delete ;
run;

proc sort data = CC nodups ;
by _all_ ;
run;

/* merge into mlst_w_pairs from identify_if_pairs_have_ST_matches.sas */
proc sort data = CC ;
by ST ;
proc sort data = mlst_w_pairs ;
by ST ;
run;

data ST_CC_relapse ;
merge CC (in=in1) mlst_w_pairs (in=in2) ;
by ST ;
if in2 ;
run;

/* no CC group for ST 840 and no CC group for ST 1181 - see documentation for picking CC */

data ST_CC_relapse2 ;
set  ST_CC_relapse ;
if ST = 840 then CC = 5 ;
if ST = 1181 then CC = 8 ;
run ;

proc freq data = ST_CC_relapse2 ;
tables CC ;
run;
    /*  
CC    Frequency
----------------
 1           4
 5          25
 8          40
15           6
30           2
45           2  */

data CC_refs ;
set ref.reference_table_locked ;
keep ST CC strain taxonomy_ID  ;
if ST = . then delete ;
if CC = . then delete ;
rename strain = reference ;
run;

data CC_refs2 ;
set CC_refs ;
if CC = 1 and reference = "MSSA476" then flag_ref = 1 ;
else if CC = 5 and reference = "ED98" then flag_ref = 1 ;
else if CC = 8 and reference = "Newman" then flag_ref = 1 ;
else if CC = 15 and reference = "ST20130941" then flag_ref = 1 ;
else if CC = 30 and reference = "TCH60" then flag_ref = 1 ;
else if CC = 45 and reference = "CA_347" then flag_ref = 1 ;
else flag_ref = 0;
run;

data CC_refs3 ;
set CC_refs2 ;
where flag_ref = 1 ;
run;

data relapse.CC_refs;
set CC_refs3 ;
run ;

proc export data = relapse.CC_refs
outfile = "!MCLAB/staph/seoung_ho_project/output/CC_refs.csv"
dmbs = csv replace ;
run;

proc export data = relapse.CC_refs
outfile = "!MCLAB/staph/seoung_ho_project/output/CC_refs_noHeader.csv"
dmbs = csv replace ;
putnames  = no ;
run;


proc sort data = CC_refs3 ;
by CC  ;
proc sort data = ST_CC_relapse2 ;
by CC ;
run ;

data ST_CC_ref_relapse ;
merge CC_refs3 (in=in1) ST_CC_relapse2 (in=in2) ;
by CC ;
if in2 ;
run ;

data relapse.sampleID_ST_CC_ref ;
retain sampleID sampleNum ST CC flag_failed_library reference taxonomy_ID;
set ST_CC_ref_relapse ;
run ;

/* export for scp to hpc for read alignments to references */
proc export data = relapse.sampleID_ST_CC_ref  
outfile = "!MCLAB/staph/seoung_ho_project/output/sampleID_ST_CC_ref.csv"
dbms = csv replace ;
run ;

proc export data = relapse.sampleID_ST_CC_ref  
outfile = "!MCLAB/staph/seoung_ho_project/output/sampleID_ST_CC_ref_noHeader.csv"
dbms = csv replace ;
putnames = no ;
run ;




