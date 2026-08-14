



%macro importing (geno) ;

proc import datafile = "!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_&geno._ujc.csv"
out = flag_5species_2_&geno. 
dbms = csv replace ;
guessingrows=max;
run ;

data flag2_5species_2_&geno. ;
set flag_5species_2_&geno.  ;
rename geneID = &geno._geneID ;
run ;

proc sort data = flag2_5species_2_&geno. ;
by &geno._jxnHash ;
run;


%mend ;

%importing (dmel6) ;
/*75986*/

add dmel6_geneid label!

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

run;

proc freq data=mel_anno;
tables flag_sexDet_gene;
run;

data mel_anno2 ;
set mel_anno ;
if flag_sexDet_gene ne 1 then flag_sexDet_gene = 0 ;
rename primary_FBgn = dmel6_geneid ;
rename symbol = mel_symbol ;
run;

proc sort data=flag2_5species_2_dmel6;
by dmel6_geneID;
proc sort data=mel_anno2;
by dmel6_geneID;

data fivespecies_dmel_with_geneid;
merge flag2_5species_2_dmel6 mel_anno2;
by dmel6_geneid;
run;

