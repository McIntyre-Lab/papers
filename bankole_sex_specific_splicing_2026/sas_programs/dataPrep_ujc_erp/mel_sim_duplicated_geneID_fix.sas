

/*

lmm has list of 'missing' dmel and dsim geneIDs that datafile have 'wrong' duplicated geneID

lmm email 12 apr 26:
I have results in the gene files with no genesetid.
I took the geneid's for the network from  the geneid's linked to components in component_geneID_edges_fiveSpecies.csv

there are 3 dmel and 16 dsim - we found in duplicated gene list so datafile gffcompare has 'wrong' geneID

output:
/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/
    dmel6_duplicated_geneID_fix.csv
    dsim2_duplicated_geneID_fix.csv
*/


proc import datafile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_dmel6_anno_files/dmel650_duplicated_genes.csv"
out = mel
dbms = csv replace ;
guessingrows = MAX ;
run;

data mel2 ;
set mel ;
keep geneID_orig geneID_new ;
run ;
proc sort data = mel2 nodups ;
by _all_ ;
run;


proc import datafile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_dsim2_anno_files/dsim202_duplicated_genes.csv"
out = sim
dbms = csv replace ;
guessingrows = MAX ;
run;

data sim2 ;
set sim ;
keep geneID_orig geneID_new ;
run ;
proc sort data = sim2 nodups ;
by _all_ ;
run;

data lmm_mel;
    length geneID_orig $15;
    input geneID_orig $;
    datalines;
FBgn0085193
FBgn0261843
FBgn0287594
;
run;

data mel_fix ;
merge lmm_mel (in=in1) mel2 (in=in2) ;
by geneID_orig ;
if in1 ;
run ;

data mel_fix2 ;
set mel_fix ;
rename geneID_orig = lmm_missing_geneID ;
run;


data lmm_sim;
    length geneID_orig $15;
    input geneID_orig $;
    datalines;
FBgn0184795
FBgn0185916
FBgn0191573
FBgn0193071
FBgn0194257
FBgn0195183
FBgn0196490
FBgn0268461
FBgn0268506
FBgn0268568
FBgn0268691
FBgn0268999
FBgn0269029
FBgn0269085
FBgn0269373
FBgn0269787
;
run ;

data sim_fix ;
merge lmm_sim (in=in1) sim2 (in=in2) ;
by geneID_orig ;
if in1 ;
run ;


data sim_fix2 ;
set sim_fix ;
rename geneID_orig = lmm_missing_geneID ;
run;

proc export data = mel_fix2 
outfile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/dmel6_duplicated_geneID_fix.csv"
dbms = csv replace ;
run;

proc export data = sim_fix2 
outfile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/dsim2_duplicated_geneID_fix.csv"
dbms = csv replace ;
run;


