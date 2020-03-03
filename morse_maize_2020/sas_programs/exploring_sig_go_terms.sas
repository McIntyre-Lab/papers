libname B73v4 "!MCLAB/useful_maize_info/RefGenV4/sasdata";

libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";



/*

GO enrichment for maize genes  -- looking at go terms across genotypes and all 5 
*/



/* create freq table of go_terms */
%macro goCnts (aspect) ;

data  &aspect ;
set tappas.go_enrichment_&aspect ;
if fisher_fdr_p le 0.01 then flag_go_sig = 1 ;
    else flag_go_sig = 0;
run ;

proc freq data = &aspect ;
where flag_go_sig = 1 ;
**by sigVarName ;
tables go_term / out = &aspect._cnts;
run;

data go_term_cnts_de_&aspect.2 ;
set &aspect._cnts ;
drop percent ;
run ;

proc sort data = go_term_cnts_de_&aspect.2 ;
by go_term ;
run;

data go_term_cnts_de_&aspect. ;
set go_term_cnts_de_&aspect.2 ;
label count =  "count_&aspect." ;
rename count = count_&aspect. ;

run ;
%mend ;

%goCnts (molFunction) ;
%goCnts (bioProcess) ;
%goCnts (cellComponent) ;

data go_term_counts ;
merge go_term_cnts_de_molFunction go_term_cnts_de_bioProcess  go_term_cnts_de_cellComponent;
by go_term ;
run;

data tappas.go_enrich_sig_go_term_table ;
set  go_term_counts ;
label  count_molFunction = "count_molFunction";
label  count_bioProcess = "count_bioProcess";
label  count_cellComponent = "count_cellComponent";
array change _numeric_ ;
do over change ;
if change =. then change =0;
end ;
run;

proc sort data = tappas.go_enrich_sig_go_term_table ;
by go_term ;
run ;

proc export data = tappas.go_enrich_sig_go_term_table
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/go_enrich_sig_go_term_table.tsv"
dbms = tab replace ;
run ;



/* create table of which genotypes have sig go terms by genotype (and all5) */
%macro tables (aspect, shrt) ;

data &aspect._go_sig ;
set &aspect. ;
where flag_go_sig = 1 ;
keep sigVarName go_id go_term flag_go_sig;
run;

proc sort data = &aspect._go_sig ;
by go_term ;
run;

proc transpose data = &aspect._go_sig out = &aspect._go_flip ;
by go_term ;
var flag_go_sig ;
id sigVarName ;
run ;

data &aspect._go_sbys2 ;
retain go_term flag_de_all5 flag_de_b73 flag_de_c123 flag_de_hp301 flag_de_mo17 flag_de_nc338 ;
set &aspect._go_flip ;
array change _numeric_ ;
do over change ;
if change =. then change =0;
end ;
drop _name_ ;
run;

data &aspect._go_sbys ;
set  &aspect._go_sbys2 ;
rename flag_de_all5 = flag_go_sig_all5 ;
rename flag_de_b73 = flag_go_sig_B73 ;
rename flag_de_c123 = flag_go_sig_C123 ;
rename flag_de_hp301 = flag_go_sig_Hp301 ;
rename flag_de_mo17 = flag_go_sig_Mo17 ;
rename flag_de_nc338 = flag_go_sig_NC338 ;
run;

proc sort data = &aspect._go_sbys ;
by decending flag_go_sig_all5 go_term ;
run ;

data tappas.go_enrich_goTerm_sig_&shrt.;
set &aspect._go_sbys ;
run ;

proc export data = tappas.go_enrich_goTerm_sig_&shrt.
outfile ="!MCLAB/maize_ozone_FINAL/2018/PacBio/go_enrichment/go_enrich_goTerm_sig_&shrt..tsv" 
dbms = tab replace ;
run;
%mend ;

%tables (molFunction, molF) ;
%tables (bioProcess, bioP) ;
%tables (cellComponent, cellC) ;

data chk_heat ;
set tappas.go_enrichment_bioProcess;
where go_term = "heat acclimation";
run;  
 



