

libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";

/*
want geneIDs ...

import file with link between dmel_jnxHash and geneID
 
import flag files for each species
    contains geneID, jnxHash, flags for if that jxnhash also in others...

merge flag file to sexsplic.cross_species_dser_2_dmel6
  on species jxnHash

then merge in dmel_geneKey
  on dmel jxnHash
 
*/

proc import datafile = "!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_dmel6_ujc_w_IR_flag.csv"
out = dmel6_geneKey
dbms = csv replace ;
guessingrows = MAX ;
run;

data dmel6_geneKey2;
set dmel6_geneKey;
rename geneID=dmel6_geneid;
keep dmel6_jxnHash geneID;
run;

proc sort data=dmel6_geneKey2;
by dmel6_jxnHash;
run;

data aa2 ;
set dmel6_genekey2 ;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39" or dmel6_geneID = "FBgn0266214" 
or dmel6_geneid = "FBgn0262425";
run;

data find ;
set dmel6_geneKey2 ;
by dmel6_jxnHash ;
retain count ;
if first.dmel6_jxnHash then count = 0 ;
count +1;
if last.dmel6_jxnHash then do ;
if count > 1 then flag_dup = 1;
else flag_dup = 0;
output ;
end ;
run ;

proc freq data = find  ;
tables flag_dup ;
run;



%macro importing (geno) ;

proc import datafile = "!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_&geno._ujc_w_IR_flag.csv"
out = flag_5species_2_&geno. 
dbms = csv replace ;
guessingrows=max;
run ;

data flag2_5species_2_&geno. ;
set flag_5species_2_&geno.  ;
rename geneID = &geno._geneID ;
drop dmel6_geneID ;
run ;

proc sort data = flag2_5species_2_&geno. ;
by &geno._jxnHash ;
run;


%mend ;

%importing (dmel6) ;
%importing (dsan1) ;
%importing (dser1) ;
%importing (dsim2) ;
%importing (dyak2) ;
/*%importing (dsimW) ; */



%macro mismatch_cnts (geno, shrt) ;

proc sort data = sexsplic.cross_species_&shrt._2_dmel6 ;
by &geno._jxnHash ;
run ;

/* merge in species geneID */
data CS_&shrt._w_gene ;
merge sexsplic.cross_species_&shrt._2_dmel6 (in=in1) flag2_5species_2_&geno. (in=in2) ;
by &geno._jxnHash ;
if in1 ;
run;

proc sort data = CS_&shrt._w_gene ;
by dmel6_jxnHash ;
run;

/* merge in dmel6 geneID */
data CS2_&shrt._w_gene_dmel6 ;
merge CS_&shrt._w_gene (in=in1) dmel6_geneKey2 (in=in2) ;
by dmel6_jxnHash ;
if in1 ;
run;


/*  number of genes matching */
data CS3_&shrt._w_gene_dmel6 ;
retain &geno._jxnHash dmel6_jxnHash flag_mel_jxnHash &geno._geneID dmel6_geneID ;
set CS2_&shrt._w_gene_dmel6 ;
if dmel6_geneID = "" then delete ;
run ;  /* NOT INCLUDING WHERE THERE IS NO DMEL6_GENEID */


/* Identify unique pairs of dser1 and dmel6 geneIDs 
        count number of &geno._geneIDs that go to */
proc sql;
    create table unique_pairs_&geno._2_dmel6 as
    select distinct &geno._geneID, dmel6_geneID
    from CS3_&shrt._w_gene_dmel6;
quit;

proc sql;
    create table &geno._cnts_combined as
    select 
        a.*,
        count(distinct b.dmel6_geneID) as num_dmel6_geneID,
        count(distinct c.&geno._geneID) as num_&geno._geneID
    from unique_pairs_&geno._2_dmel6 as a
    left join unique_pairs_&geno._2_dmel6 as b
    on a.&geno._geneID = b.&geno._geneID
    left join unique_pairs_&geno._2_dmel6 as c
    on a.dmel6_geneID = c.dmel6_geneID
    group by a.&geno._geneID, a.dmel6_geneID;
quit;

proc sort data= &geno._cnts_combined;
by &geno._geneId dmel6_geneID;
proc sort data= CS2_&shrt._w_gene_dmel6;
by &geno._geneid dmel6_geneID;
run;

/* this is now a many to many merge!!! shit */
data fivespecies_&geno._with_geneid;
merge CS2_&shrt._w_gene_dmel6 &geno._cnts_combined;
by &geno._geneid dmel6_geneID ;
run;

%mend ;

%mismatch_cnts (dsim2, dsim2) ;
%mismatch_cnts (dser1, dser) ;
%mismatch_cnts (dsan1, dsan) ;
%mismatch_cnts (dyak2, dyak) ;
/*%mismatch_cnts (dsimW, dsimW) ; 

data fiveSpecies_dsimW_with_geneID ;
set  fiveSpecies_dsimW_with_geneID ;
rename dsimW_jxnHash = dsim2_jxnHash ;
run;
*/

/* checking....
data aaa ;
set sexsplic.cross_species_dser_2_dmel6;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run;

data aa_2 ;
set CS2_dser_w_gene_dmel6  ;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run ;

data aa_3 ;
set CS3_dser_w_gene_dmel6;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run ;

data aa_4 ;
set fivespecies_dser1_with_geneid;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run ;


proc freq data = CS2_dser_w_gene_dmel6 ;
tables dser1_geneID  
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run ;

data aa_5 ;
set unique_pairs_dser_2_dmel6;
where dmel6_jxnHash = "0019ea64616021e14cfeef519ba46e9f524bae59a907c7c07fc78a371e77fa39";
run ;
*/

/* make perm */
%macro perm (geno) ;

data sexsplic.fivespecies_&geno._w_geneid_02amm ;
retain &geno._jxnHash dmel6_jxnhash &geno._geneID dmel6_geneID num_&geno._geneID num_dmel6_geneID ;
set fivespecies_&geno._with_geneid;
run ;

proc export data = sexsplic.fivespecies_&geno._w_geneid_02amm 
oufile = "!MCLAB/sex_specific_splicing/cross_species_link_files/fivespecies_&geno._w_geneid_02amm.csv"
dbms = csv replace ;
run ;

%mend ;

%perm (dser1) ;
%perm (dsan1) ;
%perm (dsim2) ;
/* %perm (dsimW) ; */
%perm (dyak2) ;



ods pdf file="!MCLAB/sex_specific_splicing/cross_species_link_files/gene_pair_jxnHash_cnts_for_geno_2_dmel6.pdf";

%macro cnts (geno) ;
title "num geneID pairs by jxnHash for &geno. and dmel6 ";
proc freq data = sexsplic.fivespecies_&geno._w_geneid_02amm noprint ;
tables num_&geno._geneID * num_dmel6_geneID / out = cnts_4_&geno. ;
run ;

proc print data = cnts_4_&geno. ; run;
title "";
%mend ;


%cnts (dser1) ;
%cnts (dsan1) ;
%cnts (dsim2) ;
/*%cnts (dsimW) ; */
%cnts (dyak2) ;

ods pdf close ;


