
libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";
/*


maize short read data

drop bad samples!  see work.list for link between old and new sample names

    B73_P4_C1_Amb = B73_Amb_4
    B73_P4_C6_Amb = B73_Amb_12
    B73_P1_C7_Ele = B73_Ele_9 

*/


proc sort data  = make.cvrg_cnts_shrtRead ;
by primary_FBgn ;
run;

/* sbys with wt_tpm */
proc transpose data = make.cvrg_cnts_shrtRead out = cvrg_shrt_flip  ;
by primary_FBgn ;
id sampleID ;
var wt_tpm ;
run;

data cvrg_shrt_sbys ;
set cvrg_shrt_flip ;
drop _name_ ;
run;

/* merge in flag_analyze tpm5*/
data flag_analyze ;
set make.onCalls_shrt_gene_tpm5 ;
keep primary_FBgn flag_analyze_: ;
run;

proc sort data = flag_analyze ;
by primary_FBgn ;
proc sort data = cvrg_shrt_sbys ;
by primary_FBgn ;
run;

data sbys_tappas_tpm ;
merge flag_analyze (in=in1) cvrg_shrt_sbys (in=in2) ;
by primary_FBgn ;
if in2 ;
run ;

/* upcase variable names */
option validvarname = upcase ;
data sbys_tappas_tpm_up ;
set sbys_tappas_tpm ;
run;
option validvarname = V7 ;

/* create associated design files:  sample condition  */
data designFile;
set make.cvrg_cnts_shrtRead ;
if trt = "Amb" then condition = "Ambient" ;
if trt = "Ele" then condition = "Ozone" ;
rename sampleID = sample ;
keep sampleID condition ;
run ;

proc sort data = designFile nodups ;
by _ALL_ ;
run;  /* 120 */

data designFile_up ;
retain newSample condition ;
set designFile ;
newSample = upcase(sample) ;
drop sample ;
rename newSample = sample ;
run ;



%macro splitting (geno) ;

/* split out genotypes into separate files */
proc sql noprint ;
select name into : varnames separated by " "
from dictionary.columns
where libname = "WORK" and memname = "SBYS_TAPPAS_TPM_UP"
    and name contains ("&geno.") or upcase(name) contains "PRIMARY_FBGN" ;
quit;

data sbys_&geno._4_tappas_tpm ;
set sbys_tappas_tpm_up (keep = &varnames) ;
where flag_analyze_&geno. = 1 ;
drop flag_analyze_&geno. ;
run ;


/* split design file into genotypes */
data df_&geno._tappas_tpm ;
set designFile_up ;
if find(sample, "&geno.") ge 1 ;
run;

proc export data = sbys_&geno._4_tappas_tpm
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/make_files_4_tappas_amm/sbys_&geno._4_tappas_tpm.tsv"
dbms = tab replace ;
run;

proc export data = df_&geno._tappas_tpm 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/make_files_4_tappas_amm/df_&geno._tappas_tpm.tsv"
dbms = tab replace ;
run;

%mend ;

%splitting (B73) ;
%splitting (MO17) ;
%splitting (C123) ;
%splitting (HP301) ;
%splitting (NC338) ;


/* comparing b73 amb to other 4 genotypes amb tpm */ 

/* merge in flag_&geno._amb_on */
data flag_amb_on ;
set make.onCalls_shrt_gene_tpm5 ;
keep primary_FBgn flag_b73_amb_on flag_mo17_amb_on flag_c123_amb_on flag_hp301_amb_on flag_nc338_amb_on ;
run;

proc sort data = flag_amb_on ;
by primary_FBgn ;
proc sort data = cvrg_shrt_sbys ;
by primary_FBgn ;
run;


data sbys_tappas_amb_tpm ;
merge flag_amb_on (in=in1) cvrg_shrt_sbys (in=in2) ;
by primary_FBgn ;
if in2 ;
run ;

/* upcase variable names */
option validvarname = upcase ;
data sbys_tappas_amb_tpm_up ;
set sbys_tappas_amb_tpm ;
run;
option validvarname = V7 ;



%macro splitting2 (geno1, geno2) ;

/* split out b73 amb and X amb genotypes into files */
proc sql noprint ;
select name into : varnames separated by " "
from dictionary.columns
where libname = "WORK" 
    and memname = "SBYS_TAPPAS_AMB_TPM_UP" 
    and (name contains ("&geno1.") or name contains ("&geno2.")) 
    and upcase(name) not contains ("ELE") 
    or upcase(name) contains "PRIMARY_FBGN"  ;
quit;


data sbys_&geno1._&geno2._amb_4_tappas_tpm ;
set sbys_tappas_amb_tpm_up  (keep = &varnames) ;
where flag_&geno1._amb_on = 1 or flag_&geno2._amb_on = 1 ;
  * drop  flag_&geno1._amb_on flag_&geno2._amb_on ;
run ;



/* split design file into genotypes */
data df_&geno1._&geno2._amb_4_tappas_tpm ;
retain sample geno ;
set designFile_up ;
if find(sample, "&geno1.") ge 1 or find(sample, "&geno2.") ge 1 ;
geno = compress(scan(sample, 1, '_')) ;
drop condition ;
rename geno = condition ;
run;

proc export data = sbys_&geno1._&geno2._amb_4_tappas_tpm
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/make_files_4_tappas_amm/sbys_&geno1._&geno2._amb_4_tappas_tpm.tsv"
dbms = tab replace ;
run;

proc export data = df_&geno1._&geno2._amb_4_tappas_tpm 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/make_files_4_tappas_amm/df_&geno1._&geno2._amb_4_tappas_tpm.tsv"
dbms = tab replace ;
run;


%mend ;

%splitting2 (B73, MO17) ;


%splitting2 (B73, C123) ;
%splitting2 (B73, HP301) ;
%splitting2 (B73, NC338) ;

