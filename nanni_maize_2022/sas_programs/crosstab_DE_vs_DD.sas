


libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

libname anno "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/* 
cross tab / table between genes differentially expressed (tappas) vs genes differentially detected (DD) 

where flags are number of genotypes each gene is DD or DE

input:
    tappas.de_genes_w_go
    pacbio.Diff_detect_amb_ele_gene 
    pacbio.pacbio.fsm_ism_isoform_zmtr

output:

*/

data DD ;
set pacbio.Diff_detect_amb_ele_gene ;
run ;

/* merge ZMgn to DD */
data wanno_4_genes;
set pacbio.fsm_ism_nic_nnc_wanno;
keep geneID ZMgn ;
run ;

proc sort data = wanno_4_genes nodups ;
by _all_ ;
proc sort data = wanno_4_genes;
by geneID ;
run ;

proc sort data = DD ;
by geneID ;
run ;

data DD_zmgn ;
merge DD (in=in1) wanno_4_genes (in=in2) ;
by geneID ;
if in1 ;
run ;

data DD_zmgn ;
retain zmgn ;
set  DD_zmgn ;
drop B73 C123 Hp301 Mo17 NC338  ;
run ;


data DE ;
format ZMgn $44.;
set tappas.de_genes_w_go ;
DE_geno_count = sum(flag_de_b73, flag_de_c123, flag_de_hp301, flag_de_mo17, flag_de_nc338) ;
ZMgn = gene;
drop go_:  gene;
run;

proc sort data = DE ;
by ZMgn ;
proc sort data = DD_zmgn ;
by ZMgn ;

data DD_and_DE ;
merge DE (in=in1) DD_zmgn (in=in2) ;
by ZMgn ;
run ;

proc sort data = DD_and_DE ;
by geneID ;
run ;

ods pdf file = "!MCLAB/maize_ozone_FINAL/pacbio_paper/penultimate_version/draft_figs_tables/DD_and_DE_crosstab.pdf" ;
proc freq data = DD_and_DE ;
tables DE_geno_count * DD_geno_count / out = DD_DE_crosstab ;
run;
ods pdf close ;




