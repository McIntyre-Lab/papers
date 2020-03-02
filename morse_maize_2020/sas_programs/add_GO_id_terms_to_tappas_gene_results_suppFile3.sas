libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname B73v4 "!MCLAB/useful_maize_info/RefGenV4/sasdata";

/*    

add go terms to suppl table 3
    currently contains:  gene, geneName, tappas results
*/

/* file3 */
data file3 ;
set tappas.tappas_results_genes_wanno ;
run;

/* goids */
proc import datafile = "!MCLAB/useful_maize_info/RefGenV4/maize.B73.AGPv4.aggregate_02amm.gaf"
out = gaf
dbms = tab replace ;
guessingrows = MAX;
run;

proc freq data = gaf ;
tables aspect ;
run;

/* split out aspects */

data go_MF go_P go_C oops ;
format term_accession $12. ;
set gaf ;
if aspect = "F" then output go_MF ;
else if aspect = "P" then output go_P ;
else if aspect = "C" then output go_C ;
else output oops ;
rename db_object_id = gene ;
rename db_object_symbol = symbol ;
run;  /* 0 in oops */

data goid_mf ;
set go_mf ;
keep gene term_accession ;
rename term_accession = GO_ID ;
run ;
proc sort data = goid_mf nodups out = dupID ;
by _all_ ;
run ;  /* no dups */


/* go terms */
proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/GO_table.txt"
out = table 
dbms = tab replace ;
guessingrows = MAX ;
run;

data goterms ;
set table ;
label go_id = "GO_ID"; 
keep GO_ID term definition ontology ;
where ontology = "MF";
rename term = GO_Term ;
run;

proc sort data = goterms nodups dupout = dups ;
by _all_ ;
run;
/* go_id not unique */

data ck ;
set table ;
where go_id = "GO:0000014";
run;  /* looks ok to dump dups */

proc freq data = goterms ;
tables go_id / out = cnts ;
run;
data check ;
set cnts ;
where count ne 1 ;
run;  /*go-id is unique */

proc sort data = goterms ;
by go_id ;
run;

/* merge goIDs with goterms */
proc sort data = goid_mf ;
by GO_ID ; 
proc sort data = goterms ;
by GO_ID ; 
run;

data go_id_term ;
merge goid_mf (in=in1) goterms (in=in2) ;
by go_id ;
if in1 then output go_id_term ;
run ;

/* merge with tappas genes */
proc sort data = go_id_term ;
by gene ;
proc sort data = gene ;
by gene ;
run ;

data add_go ;
merge gene (in=in1) go_id_term (in=in2) ;
by gene ;
if in1 then output add_go ;
run ;

data tappas.tappas_results_gene_wGO ;
retain gene geneName regulated_b73 regulated_c123 regulated_hp301 regulated_mo17 regulated_nc338 ;
set add_go ;
drop name_description ;
run ;

proc sort data = tappas.tappas_results_gene_wGO ;
by descending regulated_b73 ;
run ;

proc export data = tappas.tappas_results_gene_wGO
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/supplementary_data/Supplementary_File_3_02amm.tsv"
dbms = tab replace ;
run ;









