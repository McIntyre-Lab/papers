


libname B73v4 "!MCLAB/useful_maize_info/RefGenV4/sasdata";

libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";





/* go enrichment for genes DE in  only 1 of the 5 genotypes 
    based on code in distribution_gene_transcript_only_1_geno.sas
*/


data gene_cnts2 ;
set tappas.tappas_results_genes ;
keep gene regulated_: ;
run;

data gene_cnts ;
set gene_cnts2 ;
if regulated_B73 ne ""  then flag_on_B73 = 1 ;    else flag_on_B73 = 0 ;
if regulated_Mo17 ne "" then flag_on_Mo17 = 1 ;   else flag_on_Mo17 = 0 ;
if regulated_C123 ne "" then flag_on_C123 = 1 ;   else flag_on_C123 = 0 ;
if regulated_Hp301 ne "" then flag_on_Hp301 = 1 ; else flag_on_Hp301 = 0 ;
if regulated_NC338 ne "" then flag_on_NC338 = 1 ; else flag_on_NC338 = 0 ;
sum_on = (flag_on_B73 + flag_on_Mo17 + flag_on_C123 + flag_on_Hp301 + flag_on_NC338) ;
drop regulated_:;
run ;

proc sort data = gene_cnts ;
by sum_on ;
run;

data gene_sig_only_1_geno ;
set gene_cnts ;
if sum_on = 1 then do ;
if flag_on_B73 = 1 then flag_sig_B73only = 1 ; else flag_sig_B73only = 0 ;
if flag_on_C123 = 1 then flag_sig_C123only = 1 ; else flag_sig_C123only = 0 ;
if flag_on_Hp301 = 1 then flag_sig_Hp301only = 1 ; else flag_sig_Hp301only = 0 ;
if flag_on_Mo17 = 1 then flag_sig_Mo17only = 1 ; else flag_sig_Mo17only = 0 ;
if flag_on_NC338 = 1 then flag_sig_NC338only = 1 ; else flag_sig_NC338only = 0 ;
end ;
else if sum_on ne 1 then do ;
flag_sig_B73only = 0 ;
flag_sig_C123only = 0 ;
flag_sig_Hp301only = 0 ;
flag_sig_Mo17only = 0 ;
flag_sig_NC338only = 0 ;
end ;
drop flag_on_: sum_on  ;
run;



data GO ;
set b73v4.GO_IDs_catted ;
rename ID  = gene ;
run ;

proc sort data = gene_sig_only_1_geno ;
by gene ;
proc sort data = go;
by gene ;
run ;

data de_go ;
merge gene_sig_only_1_geno (in=in1) go (in=in2) ;
by gene ;
if in1 then output de_go ;
run ;

data DE_go_drop_empty ;
set de_go ;
if go_molFunction ne '' or go_bioprocess ne '' or go_cellcomponent ne '' ;
run;

data tappas.DE_genes_1_genotypeOnly_w_GO ;
set DE_go_drop_empty ;
run ;
/* gene set enrichment in JMP Genomics 
    
(a) 3 functional groups (= category variables,  do 1 at a time):
		(1) go_molfunction, (2) go_bioProcess and (3) go_cellComponent

(b) binary significance variables:
	flag_sig_b73only
	flag_sig_C123only
	flag_sig_Hp301only
	flag_sig_Mo17only
	flag_sig_NC338only
(c) "larger is more significant"
(d) Fisher enrichment test
(e) sig cutoff = 1 
(f) low hit threshold = 5
    high hit threshold = 100
    max category length = 512
    
data into JMP:  tappas.DE_genes_1_genotypeOnly_w_GO

data out of JMP:    (1) GO_enrich_jmp_output_1Geno_Mol
                    (2) GO_enrich_jmp_output_1Geno_Bio 
                    (3) GO_enrich_jmp_output_1Geno_Cell

*/






/* import table with GO names */

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/GO_table.txt"
out = table 
dbms = tab replace ;
guessingrows = MAX ;
run;

data table2 ;
set table ;
keep go_id term  ontology ;
run;

proc sort data = table2 nodups  ;
by _all_ ;
run;
/* go_id not unique */

proc freq data = table2 ;
tables go_id / out = cnts ;
run;
data check ;
set cnts ;
where count ne 1 ;
run;  /*go-id is unique */

proc sort data = table2 ;
by go_id ;
run;



%macro goStuff (aspect) ;

data go_&aspect ;
retain sigVarName BsigVarName go_id ;
set tappas.GO_enrich_jmp_output_1Geno_&aspect ;
where c0_less_than_c1 = 1;
go_id = tranwrd(category, "go:", "GO:") ;
keep   sigVarName
 Fisher_FDR_p
 Fisher_Raw_p
 go_id
 _xval_
  _x_over_m_
_mval_
 _yval_
 _y_over_n_
 _nval_ ;
 run; 

proc sort data = go_&aspect ;
by go_id;
run;

data go_id_&aspect ;
merge go_&aspect (in=in1) table2 (in=in2) ;
by go_id ;
if in1 ;
run;

proc sort data = go_id_&aspect ;
by Fisher_FDR_p ;
run ;

%mend ;

%goStuff (mol) ;
%goStuff (bio) ;
%goStuff (cell) ;

%macro fixup (short, long) ;

data tappas.go_enrich_1Geno_&long. ;
retain sigVarName go_id  term ontology ;
set go_id_&short ;
label go_id = "GO_id" ;
rename term = GO_term ;
rename ontology = GO_ontology ;
run;

proc sort data = tappas.go_enrich_1Geno_&long;
by fisher_fdr_p ;
run ;

proc export data = tappas.go_enrich_1Geno_&long. 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/go_enrichment_1Geno_DE_genes_&long..csv"
label 
dbms = csv replace ;
run ;

%mend ;

%fixup (mol, molFunction) ;
%fixup (bio, bioProcess) ;
%fixup (cell, cellComponent) ;







