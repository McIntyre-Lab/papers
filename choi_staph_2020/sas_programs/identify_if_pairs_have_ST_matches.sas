

libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
filename mymacros "!MCLAB/maize_ozone_final/2014/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);


/*
do isolate pairs have same stringMLST ST data?

import stringMLST results
merge in sampleID

*/

proc import datafile = "/home/ammorse/TB14/staph_relapse/stringMLST_analysis/stringMLST_out.tsv"
out = mlst
dbms = tab replace ;
guessingrows = MAX ;
run;

data mlst2 ;
retain sampleNum ;
set mlst ;
sampleNum = scan(sample, 2, '-');
run ;

proc sort data = mlst2 ;
by sampleNum ;
run;
proc sort data = relapse.sampleNum_sampleID ;
by sampleNum ;
run;

data mlst3 ;
merge mlst2 (in=in1) relapse.sampleNum_sampleID (in=in2) ;
by sampleNum ;
run ;

data mlst4 ;
retain sampleNum sampleID ST;
set mlst3 ;
keep sampleID ST sampleNum flag_failed_library ;
run;

data pairs ;
set relapse.isolate_pairs ;
pairNum = _N_;
keep isolate_1 isolate_2 pairNum ;
run ;

data pair1 ;
set pairs ;
keep isolate_1 pairNum ;
rename isolate_1 = sampleID;
run;

data pair2 ;
set pairs ;
keep isolate_2 pairNum ;
rename isolate_2 = sampleID;
run;

data pairs_stack ;
set pair1 pair2 ;
run ;

proc sort data = pairs_stack ;
by sampleID  ;
proc sort data = mlst4 ;
by sampleID ;
run;

data mlst_w_pairs ;
merge mlst4 (in=in1) pairs_stack (in=in2) ;
by sampleID ;
run ;

proc sort data = mlst_w_pairs ;
by pairNum ;
run ;


proc transpose data = mlst_w_pairs out = flip1 prefix = ST_isolate_;
by pairNum ;
var ST ;
run ;

data ST_for_relapse_pairs ;
set flip1 ;
if ST_isolate_1 = ST_isolate_2 then flag_pair_ST_match = 1 ;
else if ST_isolate_1 ne ST_isolate_2 then flag_pair_ST_match = 0 ;
else flag_pair_ST_match = 2 ;
run;

proc freq data = ST_for_relapse_pairs ;
tables flag_pair_ST_match ;
run;
    /* 1 mismatch - this is one where library failed */

    data find ;
    set ST_for_relapse_pairs ;
    where flag_pair_ST_match ne 1 ;
    run;

/* merge in sampleIDs */
proc sort data = pairs ;
by pairNum ;
proc sort data = ST_for_relapse_pairs ;
by pairNum ;
run ;

data relapse.ST_for_relapse_pairs  ;
merge pairs (in = in1) ST_for_relapse_pairs (in=in2) ;
by pairNum ;
drop _name_ ;
run ;

proc export data = relapse.ST_for_relapse_pairs 
outfile = "!MCLAB/staph/seoung_ho_project/output/ST_for_relapse_pairs.tsv"
dbms = tab replace ;
run;









