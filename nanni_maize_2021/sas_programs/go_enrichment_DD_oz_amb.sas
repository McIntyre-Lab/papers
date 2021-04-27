
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

libname anno "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/*

go enrichment of differentially detected genes 

    detected in x condition in 3 or more genotypes
    merge in goIDs
    export for jmp
    
*/




%macro DD_Prep (trt) ;

data dd_&trt._gene_4_go ;
format geneID $44. ;
set pacbio.diff_detect_&trt._gene ;
where number_genotypes_detected ge 3;
if B73 ne "" then flag_DD_&trt._B73 = 1 ; else flag_DD_&trt._B73 = 0 ;
if C123 ne "" then flag_DD_&trt._C123 = 1 ; else flag_DD_&trt._C123 = 0 ;
if Hp301 ne "" then flag_DD_&trt._Hp301 = 1 ; else flag_DD_&trt._Hp301 = 0 ;
if Mo17 ne "" then flag_DD_&trt._Mo17 = 1 ; else flag_DD_&trt._Mo17 = 0 ;
if NC338 ne "" then flag_DD_&trt._NC338 = 1 ; else flag_DD_&trt._NC338 = 0 ; 
keep geneID flag_DD_: ;
run ;

proc sort data = dd_&trt._gene_4_go ;
by geneID ;
run;

%mend ;

%DD_Prep (ele) ;
%DD_Prep (amb) ;

/* on calls where at least 1 is on */
proc contents data = pacbio.sub_geno_trt_gene_onCall_tpm0;
run ;

data onCalls ;
format geneID $ 44. ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
keep geneID ;
run ;

proc sort data = onCalls ;
by geneID ;
run;

data for_GO ;
merge onCalls dd_ele_gene_4_go dd_amb_gene_4_go ;
by geneID ;
run ;

/* check how many drop if don't include ones off in both */
data oncalls_chk ;
set pacbio.sub_geno_trt_gene_onCall_tpm0;
if flag_B73_trt_trnscpt_on0 = 0 and flag_c123_trt_trnscpt_on0 = 0 and flag_Hp301_trt_trnscpt_on0 = 0 and flag_Mo17_trt_trnscpt_on0 = 0 and flag_NC338_trt_trnscpt_on0 = 0 ;
run;  /* only 22 where off in both conditions for all genotypes */


/* set all missing values to 0 */
data DD_for_GO ;
set for_GO ;
array change _numeric_ ;
    do over change;
    if change = . then change = 0;
    end ;
run;


/* check make sure 1's all good */
title "after" ;
proc freq data = DD_for_GO ;
tables flag_DD_: ;
run;
title "before" ;
proc freq data = for_GO ;
tables flag_DD_:  ;
run;
title "";

/* merge in ZMgn  */

data wanno_4_genes;
set pacbio.pacbio.fsm_ism_isoform_zmtr;
drop isoform ;
run ;

proc sort data = wanno_4_genes nodups ;
by _all_ ;
proc sort data = wanno_4_genes;
by geneID ;
run ;

proc sort data = DD_for_GO ;
by geneID ;
run ;

data DD_ge_3geno ;
merge DD_for_GO (in=in1)  wanno_4_genes ;
by geneID ;
if in1 ;
run;

data DD_ge_3geno ;
retain geneID zmgn geneName ;
set DD_ge_3geno ;
run;
/*
proc freq data = DD_ge_3geno ;
tables geneID / out = cnts ;
run;
data ck_cnts ;
set cnts ;
where count ne 1 ;
run;

data flags ;
set ck_cnts ;
flag_mult = 1;
keep geneID flag_mult ;
run;

data DD_3geno ;
merge DD_ge_3geno flags ;
by geneID ;
run;

data DD2_3geno ;
set DD_3geno ;
where flag_mult = 1 ;
run;
*/


/* merge in go IDs */
data GO ;
set anno.GO_IDs_catted ;
rename ID  = zmgn ;
run ;

proc sort data = DD_ge_3geno ;
by zmgn ;
proc sort data = go;
by zmgn ;
run ;

data dd_go ;
merge DD_ge_3geno (in=in1) go (in=in2) ;
by zmgn ;
if in1 then output dd_go ;
run ;

data DD_ge_3geno_4_go ;
set dd_go ;
if go_molFunction ne '' or go_bioprocess ne '' or go_cellcomponent ne '' ;
drop geneID geneName ;
run ;

data pacbio.DD_genes_w_GOIDs ;
set DD_ge_3geno_4_go ;
run ;

/*gene set enrichment in JMP Genomics 
    
(a) 3 functional groups (= category variables):
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
set pacbio.go_enrich_jmp_output_DD_&aspect ;
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
%goStuff (cel) ;


%macro fixup (short, long) ;

data pacbio.go_enrich_DD_ozamb_&long. ;
retain sigVarName go_id  term ontology ;
set go_id_&short ;
label go_id = "GO_id" ;
rename term = GO_term ;
rename ontology = GO_ontology ;
run;

proc sort data = pacbio.go_enrich_DD_ozamb_&long.;
by sigVarName fisher_fdr_p ;
run ;

proc export data = pacbio.go_enrich_DD_ozamb_&long. 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/go_enrich_DD_oz_amb_&long..csv"
label 
dbms = csv replace ;
run ;

%mend ;

%fixup (mol, molFunction) ;
%fixup (bio, bioProcess) ;
%fixup (cel, cellComponent) ;











