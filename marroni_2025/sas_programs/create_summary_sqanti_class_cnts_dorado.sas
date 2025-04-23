

filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);

/*

create single classifcaion count file

*/


/* import design file with sampleID */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv"
out = df
dbms = csv replace ;
guessingrows = max ;
run ;

data design ;
retain sample sample_w_HL sampleID TR tech ;
length tech $6. ;
set df ;
TR = compress('TR'||runID) ;
sampleID = compress(sample||'_TR'||runID) ;
*keep sampleID ;
if find(runName, "prom") ge 1 then tech = "prom" ;
else if find(runName, "first") ge 1 then tech = "minion" ;
else if find(ONT_dir, "McIntyre") ge 1 then tech = "prom" ;
else if find(ONT_dir, "_FAR") ge 1 then tech = "minion" ;
else tech = "oops" ;
run ;

data design_sampleID ;
set design ;
keep sampleID ;
run;

proc sort data = design_sampleID nodups ;
by _all_ ;
run;


/* import sqanti class counts */
%macro imp_out (type) ;

%macro imp (sampleID) ;

proc import datafile = "/TB14/TB14/pxs_lmm_dros_data/dorado_basecalling/sqanti_cnts/&sampleID._&type._class_cnt.csv"
out = c_&sampleID._&type.
dbms = csv replace ;
run ;

data c2_&sampleID._&type. ;
set c_&sampleID._&type.  ;
sampleID = "&sampleID" ;
type = "&type";
run;

proc transpose data = c2_&sampleID._&type. out = cc_&sampleID._&type. prefix = &type._ ;
by sampleID type ;
id structural_category;
var count ;
run;

data &type._&sampleID._cc ;
set cc_&sampleID._&type. ;
drop type _name_ ;
run ;

%mend ;

%iterdataset(dataset=design_sampleID, function=%nrstr(%imp(&sampleID);)); 
/* %imp  (dm11037_01h_rep1_TR2); */

%mend ;
%imp_out (fl_rsc) ;
%imp_out (unclass) ;


data all_fl_rsc ;
set fl_rsc_: ;
run;

data all_unclass ;
set unclass_: ;
run;


data sqanti_class_cnts ;
merge all_: ;
by sampleID ;
run ;

/*  sum fl_rsc_FSM + unclass_FSM 
    is raw_flip_FSM > the above sum?  */

/* below for ISM NNC NIC */
data check ;
set sqanti_class_cnts ;
sum_fl_rsc_unclass_FSM = (fl_rsc_full_splice_match + unclass_full_splice_match) ;
run;


proc export data = sqanti_class_cnts 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/dorado_basecalling/sqanti_class_cnt_table_dorado.csv"
dbms = csv replace ;
run ;


