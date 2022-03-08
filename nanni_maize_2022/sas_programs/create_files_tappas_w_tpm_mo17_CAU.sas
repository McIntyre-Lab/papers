
/*


maize short read data

drop bad samples!  see work.list for link between old and new sample names

    B73_P4_C1_Amb = B73_Amb_4
    B73_P4_C6_Amb = B73_Amb_12
    B73_P1_C7_Ele = B73_Ele_9 

*/

/*
proc contents data = cvrg_cnts_shrtRead ;
run;
*/

proc sort data  = cvrg_cnts_shrtRead ;
by gene_id ;
run;

/* sbys with wt_tpm */
proc transpose data = cvrg_cnts_shrtRead out = cvrg_shrt_flip  ;
by gene_id ;
id sampleID ;
var wt_tpm ;
run;

data cvrg_shrt_sbys ;
set cvrg_shrt_flip ;
drop _name_ ;
run;

/* merge in flag_analyze tpm5*/
data flag_analyze ;
set onCalls_shrt_gene_tpm5 ;
keep gene_id flag_analyze_: ;
run;

proc sort data = flag_analyze ;
by gene_id ;
proc sort data = cvrg_shrt_sbys ;
by gene_id ;
run;

data sbys_tappas_tpm ;
merge flag_analyze (in=in1) cvrg_shrt_sbys (in=in2) ;
by gene_id ;
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
set cvrg_cnts_shrtRead ;
if trt = "Amb" then condition = "Ambient" ;
if trt = "Ele" then condition = "Ozone" ;
rename sampleID = sample ;
keep sampleID condition ;
run ;

proc sort data = designFile nodups ;
by _ALL_ ;
run;  /* 117 */

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
    and name contains ("&geno.") or upcase(name) contains "GENE_ID" ;
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
outfile = "/nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/cvr_cnts_gene_mo17_cau/sbys_&geno._4_tappas_tpm.tsv"
dbms = tab replace ;
run;

proc export data = df_&geno._tappas_tpm 
outfile = "/nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/cvr_cnts_gene_mo17_cau/df_&geno._tappas_tpm.tsv"
dbms = tab replace ;
run;

%mend ;

%splitting (B73) ;
%splitting (MO17) ;
%splitting (C123) ;
%splitting (HP301) ;
%splitting (NC338) ;

