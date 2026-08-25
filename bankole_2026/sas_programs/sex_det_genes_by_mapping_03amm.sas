

libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";



/* check for sex det genes 
    import mel sex det gene list */
proc import datafile = "!MCLAB/useful_dmel_data/gene_lists/pathways/sex_det.csv"
out = mel_sex 
dbms = csv replace ;
guessingrows = MAX ;
run;

data mel_sex2 ;
set  mel_sex ;
flag_sexDet_gene = 1 ;
run ;

proc import datafile = "!MCLAB/useful_dmel_data/flybase650/dmel_annotation/fbgn_annotation_ID.csv"
out = mel_symbol
dbms = csv replace ;
guessingrows = MAX ;
run;

data mel_symbol2 ;
set  mel_symbol ;
keep symbol primary_fbgn ;
run ;

proc sort data = mel_symbol2 nodups ;
by _all_ ;
run;

proc sort data = mel_symbol2 ;
by primary_fbgn ;
proc sort data = mel_sex2 ;
by primary_fbgn ;
run;

data mel_anno  oops;
merge mel_symbol2 (in=in1) mel_sex2 (in=in2) ;
by primary_fbgn ;
if in1 then output mel_anno ;
if flag_sexDet_gene ne 1 then flag_sexDet_gene = 0 ;
run;

data mel_anno2 ;
set mel_anno ;
rename primary_FBgn = dmel6_geneid ;
rename symbol = mel_symbol ;
run;



%macro others (geno, geno2) ;

/* check uniq */
proc sort data= sexsplic.fivespecies_&geno._w_geneid_02amm out=check nodupkey;
by &geno._geneid &geno2._jxnHash;
run; /*no duplicates*/
 
proc sort data = sexsplic.fivespecies_&geno._w_geneid_02amm ;
by dmel6_geneID ;
proc sort data= mel_anno2 ;
by dmel6_geneid;
run;

data fivespecies_&geno._geneid_sexFlag;
merge sexsplic.fivespecies_&geno._w_geneid_02amm (in=in1) mel_anno2 (in=in2);
by dmel6_geneId;
if in1;
if flag_sexDet_gene ne 1 then flag_sexDet_gene = 0;
run;

/* make perm */
data sexsplic.fivespecies_&geno._geneid_sexFlag ;
set fivespecies_&geno._geneid_sexFlag ;
run ;

proc export data = sexsplic.fivespecies_&geno._geneid_sexFlag 
outfile = "!MCLAB/sex_specific_splicing/cross_species_link_files/fivespecies_&geno._geneid_sexFlag.csv"
dbms = csv replace ;
run;

%mend ;

%others (dser1, dser1) ;
%others (dsan1, dsan1) ;
%others (dsim2, dsim2) ;
/*%others (dsimw, dsim2) ; */
%others (dyak2, dyak2) ;



ods pdf file="!MCLAB/sex_specific_splicing/cross_species_link_files/sex_det_gene_counts.pdf";

%macro counting (geno) ;

title "&geno.:  Freq flag_sexDet_gene " ;
proc freq data=sexsplic.fivespecies_&geno._geneid_sexFlag;
tables flag_sexdet_gene;
run;
title "";

title "&geno.:  Freq num_dmel6_geneID where flag_sexdet_gene = 1";
proc freq data=sexsplic.fivespecies_&geno._geneid_sexFlag;
where flag_sexdet_gene=1;
tables num_dmel6_geneid;
run;
title "";

data list2_&geno. ;
set sexsplic.fivespecies_&geno._geneid_sexFlag ;
where flag_sexDet_gene = 1 ;
run ;

proc sort data = list2_&geno. ;
by &geno._geneid;
run ;

data list_&geno. ;
set  list2_&geno. ;
by &geno._geneid;
if first.&geno._geneid;
keep &geno._geneid;
run;

proc sort data = list_&geno.  ;
by &geno._geneID ;
proc sort data = sexsplic.fivespecies_&geno._geneid_sexFlag ;
by &geno._geneID ;
run ;

data &geno._sexdet_jxns;
merge sexsplic.fivespecies_&geno._geneid_sexFlag  (in=in1) list_&geno.(in=in2);
by &geno._geneID;
if in2;
run;

data &geno._sexdet_jxns2 ;
set &geno._sexdet_jxns ;
keep mel_symbol flag_sexDet_gene num_: ;
run ;

proc sort data = &geno._sexdet_jxns2 nodups ;
by _all_ ;
run ;

title "&geno.:  sexDet mel_symbols";
proc print data = &geno._sexdet_jxns2 ; run;
title "";

%mend ;

%counting (dser1) ;
%counting (dsan1) ;
%counting (dsim2) ;
/*%counting (dsimw) ; */
%counting (dyak2) ;

ods pdf close ;






