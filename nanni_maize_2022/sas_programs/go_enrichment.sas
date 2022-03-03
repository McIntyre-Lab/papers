libname B73v4 "!MCLAB/useful_maize_info/RefGenV4/sasdata";

libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";



/*

GO enrichment for maize genes 

after GO enrichment in JMP - back into program for further analysis

*/


data DE_prep ;
set tappas.tappas_results_genes ;
if dea_result_b73 = "DE" then flag_de_b73 = 1 ;
    else flag_de_b73 = 0;
if dea_result_c123 = "DE" then flag_de_c123 = 1 ;
    else flag_de_c123 = 0;
if dea_result_hp301 = "DE" then flag_de_hp301 = 1 ;
    else flag_de_hp301 = 0;
if dea_result_mo17 = "DE" then flag_de_mo17 = 1 ;
    else flag_de_mo17 = 0;
if dea_result_nc338 = "DE" then flag_de_nc338 = 1 ;
    else flag_de_nc338 = 0;
run;

data DE_gene ;
set de_prep ;
keep gene flag_: ;
run;

data GO ;
set b73v4.GO_IDs_catted ;
rename ID  = gene ;
run ;

proc sort data = DE_gene ;
by gene ;
proc sort data = go;
by gene ;
run ;

data de_go ;
merge DE_gene (in=in1) go (in=in2) ;
by gene ;
if in1 then output de_go ;
run ;

data DE_go_drop_empty ;
set de_go ;
if go_molFunction ne '' or go_bioprocess ne '' or go_cellcomponent ne '' ;
sum_genotype_DE = (flag_de_b73 + flag_de_c123 + flag_de_hp301 + flag_de_mo17 + flag_de_nc338) ;
if sum_genotype_DE = 5 then flag_de_all5 = 1;
    else flag_de_all5 = 0 ;
run;

data tappas.DE_genes_w_GO ;
retain gene flag_de_all5   ;
set de_go_drop_empty ;
drop sum_genotype_DE ;
run ;

/* gene set enrichment in JMP Genomics 
    
(a) 3 functional groups (= category variables,  do 1 at a time):
		(1) go_molfunction, (2) go_bioProcess and (3) go_cellComponent

(b) binary significance variables:
	flag_de_all5
    flag_de_b73
    flag_de_c123
    flag_de_hp301
    flag_de_mo17
    flag_de_nc338
(c) "larger is more significant"
(d) Fisher enrichment test
(e) sig cutoff = 1 
(f) low hit threshold = 5
    high hit threshold = 100
    max category length = 512
    
data into JMP:  tappas.DE_genes_w_GO

data out of JMP:    (1) GO_enrich_jmp_output_Mol
                    (2) GO_enrich_jmp_output_Bio 
                    (3) GO_enrich_jmp_output_Cell

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


/* import gene set enrichment results */
%macro goStuff (aspect) ;

data go_&aspect ;
retain sigVarName BsigVarName go_id ;
set tappas.go_enrich_jmp_output_&aspect ;
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

data tappas.go_enrichment_&long. ;
retain sigVarName go_id  term ontology ;
set go_id_&short ;
label go_id = "GO_id" ;
rename term = GO_term ;
rename ontology = GO_ontology ;
run;

proc sort data = tappas.go_enrichment_&long;
by sigVarName fisher_fdr_p ;
run ;

proc export data = tappas.go_enrichment_&long. 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/go_enrichment_DE_genes_&long..csv"
label 
dbms = csv replace ;
run ;

%mend ;

%fixup (mol, molFunction) ;
%fixup (bio, bioProcess) ;
%fixup (cell, cellComponent) ;















