


libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";


/*

to each cross species file, add flags 1) if species_jxnHash in species data and 2) if dmel6_jxnHash in dmel6 data 
*/

%macro imp_anno (species) ;

/* import dmel6 full anno with data flag*/
proc import datafile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_&species._full_annotation_w_dataFlag.csv"
out = &species._anno 
dbms = csv replace ;
guessingrows = MAX ;
run;

data &species._anno2 ;
set &species._anno ;
keep &species._jxnHash flag_jxnHash_in_data ;
rename flag_jxnHash_in_data  = flag_&species._jxnHash_in_data ;
run ;

proc sort data = &species._anno2 nodups ;
by _all_ ;
run;

%mend ;

%imp_anno (dmel6) ;
%imp_anno (dsim2) ;
%imp_anno (dser1) ;
%imp_anno (dsan1) ;
%imp_anno (dyak2) ;


/* drop cat_dmel6_transcriptID from dmel6_anno2  */
data dmel6_anno3 ;
set dmel6_anno2 ;
*drop cat_dmel6_transcriptID ;
run;


%macro imp_cross (species) ;

/* import dsim2 cross species file */
data cross_&species. ;
set sexsplic.fivespecies_&species._w_geneid_02amm ;
*drop dmel6_geneID ;
run ;

proc sort data = cross_&species. ;
by &species._jxnHash ;
run ;

/* merge in flag_dsim2_jxnHash_in_data */
data cross_&species._f1 oops ;
merge cross_&species. (in=in1) &species._anno2 (in=in2) ;
by &species._jxnHash ;
if in1 then output cross_&species._f1 ;
else output oops ;
run;

proc sort data = cross_&species._f1 ;
by dmel6_jxnHash ;
proc sort data = dmel6_anno3 ;
by dmel6_jxnHash ;
run;

/* merge in flag_dmel6_jxnHash_in_data */
data cross_&species._f2 oops2 ;
merge cross_&species._f1 (in=in1) dmel6_anno3 (in=in2) ;
by dmel6_jxnHash ;
if in1 then output cross_&species._f2 ;
else output oops2 ;
run;

data fiveSpec_&species._geneid_dataFlags ;
retain &species._jxnHash dmel6_jxnHash flag_&species._jxnHash_in_data flag_dmel6_jxnHash_in_data cat_dmel6_transcriptID ;
set cross_&species._f2 ;
run;

%mend ;

%imp_cross (dsim2) ;
%imp_cross (dser1) ;
%imp_cross (dsan1) ;
%imp_cross (dyak2) ;


/* make perm and export */
%macro perm (species) ;

data sexsplic.fiveSpec_&species._geneid_dataFlags  ;
set fiveSpec_&species._geneid_dataFlags  ;
run ;

proc export data = fiveSpec_&species._geneid_dataFlags  
outfile = "!MCLAB//sex_specific_splicing/cross_species_link_files/fivespecies_&species._w_geneid_dataFlag.csv"
dbms = csv replace ;
run;

%mend ;

%perm (dsim2) ;
%perm (dser1) ;
%perm (dsan1) ;
%perm (dyak2) ;

/*

fivespecies_w_geneid_dataFlag file is:

FBgn0262425, which has 3 ER.

Looking at infoERP and the ER GTF, that jxnHash should have:

FBgn0266214, which correctly has 5 ER.
*/






