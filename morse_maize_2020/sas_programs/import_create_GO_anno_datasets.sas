libname B73v4 "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/*
import gaf containing GO from maize-GAMER 
*** ID and symbol are identical!

    aspects:  F = molecular function, P = biological process, C = cellular component
prep for GO enrichment
*/


proc import datafile = "!MCLAB/useful_maize_info/RefGenV4/maize.B73.AGPv4.aggregate_02amm.gaf"
out = go 
dbms = tab replace ;
guessingrows = MAX;
run;

proc freq data = go ;
tables aspect ;
run;

/* split out aspects */

data go_F go_P go_C oops ;
set go ;
if aspect = "F" then output go_F ;
else if aspect = "P" then output go_P ;
else if aspect = "C" then output go_C ;
else output oops ;
rename db_object_id = ID ;
rename db_object_symbol = symbol ;
run;  /* 0 in oops */

/* cat together go terms - want unique gene / transcript ids */
/* (1) F */
proc transpose data = go_F out = go_f_sbys ;
by ID ;
var term_accession ;
run;  /* max number is 18 */

options missing = '';
data B73v4.go_cat_molFunction ;
set go_f_sby ;
length GO_molFunction $256. ;
GO_molFunction = catx('|', of col1-col18);
drop col1-col18 _name_;
run;
options missing =.;

/* (2) P */
proc transpose data = go_P out = go_P_sbys ;
by ID ;
var term_accession ;
run;  /* max number to cat is 64 -- too long for buffer */

options missing = '';
data B73v4.go_cat_bioProcess ;
set go_P_sbys ;
length GO_bioProcess $750. ;
GO_bioProcess = catx('|', of col1-col64);
drop col1-col64 _name_ ;
run;
options missing =.;

/* (3) C */
proc transpose data = go_C out = go_C_sbys ;
by ID ;
var term_accession ;
run;  /* max number to cat is 22 -- too long for buffer */

options missing = '';
data B73v4.go_cat_cellComponent ;
set go_C_sbys ;
length GO_cellComponent $256. ;
GO_cellComponent = catx('|', of col1-col22);
drop col1-col22 _name_ ;
run;
options missing =.;

data all  ;
merge B73v4.go_cat_molFunction (in=in1) B73v4.go_cat_bioProcess (in=in2) B73v4.go_cat_cellComponent (in=in3);
by ID ;
run ;

data B73v4.GO_IDs_catted ;
set all ;
run;


