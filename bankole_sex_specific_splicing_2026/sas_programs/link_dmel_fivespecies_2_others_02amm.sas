libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";


/*get the mel gene symbols and sexdet info*/
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


proc import datafile = "!MCLAB/useful_dmel_data/gene_lists/pathways/sex_det.csv"
out = mel_sex 
dbms = csv replace ;
guessingrows = MAX ;
run;

data mel_sex2 ;
set  mel_sex ;
flag_sexDet_gene = 1 ;
run ;

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
rename primary_FBgn = geneid ;
rename symbol = mel_symbol ;
run;



/*  import fiveSpecies on dmel coord -->  link to  jxnHash for the outher four species
    import link files for all species on the dmel genome 
          -- example for dser1_2_dser1_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv:  
                has variables dser1_jxnHash and dmel6_jxnHash

    merge 5species on dmel to sim_2_dmel link file on dmel6_jxnHash
  
    for dsmel6_jxnHashes with NO link to other species_jxnHash should be the 0 in the flag file...  
    
    create flag for all dsim2_jxnHashs with link to dmel6_jxnHash

                                               Cumulative
 flag_mel_jxnHash    Frequency     Percent     Frequency
 -------------------------------------------------------
                0       11558       14.85         11558
                1       66250       85.15         77808                                      
*/

PROC IMPORT OUT= WORK.fivespecies_2dmel
            DATAFILE= "!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_dmel6_ujc_w_IR_flag.csv"
			DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;


proc sort data=fivespecies_2dmel ;
by geneid;
run;

data fivespecies_2dmel_id2;
merge fivespecies_2dmel (in=in1) mel_anno2;
by geneID;
if in1;
run;

proc sort data=fivespecies_2dmel_id2;
by dmel6_jxnhash;
run;

/*76877*/    

/* import files linking each species to dmel genome */
%import_links2 (dsan11_2_dsan1_ujc,  dmel6,   dsan1 );

%macro import_links2 (annotation, genome,  anno  );

PROC IMPORT OUT= WORK.&anno._2_&genome
            DATAFILE= "!MCLAB/sex_specific_splicing/cross_species_link_files/&annotation._2_&genome._noGeneID_ujc_xscript_link.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;


data &anno._2_&genome._id;
set &anno._2_&genome;
rename transcriptid=&anno._jxnHash
jxnHash=&genome._jxnHash;
drop geneid ; /*could keep jxnString if it helps Anna*/
run;

proc sort data=&anno._2_&genome._id;
by dmel6_jxnHash;
run;

proc transpose data=&anno._2_&genome._id out=dmel6_with_&anno (drop=_name_) prefix=&anno._;
var &anno._jxnHash;
by dmel6_jxnHash;
run;

/* playing */

proc sort data=dsan1_2_dmel6_id;
by dmel6_jxnHash dsan1_jxnHash;
run;
/*
data first ;
set dsan1_2_dmel6_id ;
by dmel6_jxnHash ;
if first.dmel6_jxnHash ;
run;
*/

proc sort data = dsan1_2_dmel6_id nodups ;
by dsan1_jxnHash ;  
run ;  /* should be uniq */

proc import datafile = "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_dsan1_anno_files/dsan11_2_dsan1_ujc_xscript_link.csv"
out = dsan11_2_dsan1_link 
dbms = csv replace ;
guessingrows = MAX ;
run;

data dsan11_2_dsan1_link2 ;
set dsan11_2_dsan1_link ;
rename jxnHash = dsan1_jxnHash ;
rename geneID = dsan1_geneID ;
drop transcriptID ;
run ;

proc sort data = dsan11_2_dsan1_link2 nodups ;
by _all_ ;
run;

proc sort data = dsan11_2_dsan1_link2 ;
by dsan1_jxnHash ;
run;

data trying_merge oops ;
merge dsan1_2_dmel6_id (in=in1) dsan11_2_dsan1_link2 (in=in2) ;
by dsan1_jxnHash ;
if in1 then output trying_merge ;
else output oops ;
run;

/* output trying merge 
create new var cat of dsan_jxnHash_dsan_geneID
nodups

do for all species 

merge each into fivespecies dmel6 full anno --> create new 5speces dmel6 with other species

check if flag_dser11_2_dser1_ujc = 1 then should be dser1_geneID

*/



data fiveSpecies_dmel6_with_dsan1;
merge fivespecies_2dmel_id2 (in=in1) first(in=in2);
by dmel6_jxnHash;
if in1 then output fiveSpecies_dmel6_with_dsan1; 
rename geneid=dmel6_geneid;
keep geneid dmel6_jxnHash dsan1_jxnHash flag_dsan11_2_dsan1_ujc ;
run;

data logic_ck ;
set fiveSpecies_dmel6_with_dsan1 ;
if flag_dsan11_2_dsan1_ujc = 1 and dsan1_jxnHash = "" then output  ;
run;  /* 0 obs */

/* */

/*
data sexsplic.fiveSpecies_dmel6_with_&anno;
merge fivespecies_2dmel_id2 (in=in1) dmel6_with_&anno(in=in2);
by dmel6_jxnHash;
if in1 then output sexsplic.fiveSpecies_dmel6_with_&anno; 
rename geneid=dmel6_geneid;
keep geneid dmel6_jxnHash flag_dmel650_2_dmel6_ujc flag_&annotation mel_symbol flag_sexDet_gene ;
run;

proc export data = sexsplic.fiveSpecies_dmel6_with_&anno
    outfile = "!MCLAB/sex_specific_splicing/cross_species_link_files/fiveSpecies_dmel6_with_&anno..csv" 
    dbms = csv replace ;
run;
*/

%mend ;

%import_links2 (dsan11_2_dsan1_ujc,  dmel6,   dsan1 );
%import_links2 (dser11_2_dser1_ujc,  dmel6,   dser1 );
%import_links2 (dyak21_2_dyak2_ujc,  dmel6,   dyak2 );

%import_links2 (dyak21_2_dyak2_ujc,  dmel6,   dyak2 );



/*sim2 and simw need a different approach*/

*%import_links2 (dsim202_2_dsim2_ujc, dmel6,   dsim2 ); /*25718 there are TWO (or more)  sim jxnhash's that will map to the same mel jxnHash;*/
*%import_links2 (dsimWXD_2_dsim2_ujc, dmel6,   dsimW ); 


PROC IMPORT OUT= WORK.dsim202_2_dsim2_ujc_2_dmel6
            DATAFILE= "!MCLAB/sex_specific_splicing/cross_species_link_files/dsim202_2_dsim2_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;

PROC IMPORT OUT= WORK.dsimw_2_dsim2_ujc_2_dmel6
            DATAFILE= "!MCLAB/sex_specific_splicing/cross_species_link_files/dsimWXD_2_dsim2_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;

data dsim_2_dmel6;
set dsim202_2_dsim2_ujc_2_dmel6 dsimw_2_dsim2_ujc_2_dmel6 ;
rename transcriptid=dsim2_jxnHash
jxnHash=dmel6_jxnHash;
drop geneid jxnString; 
run;

/*get rid of duplicates*/
proc sort data=dsim_2_dmel6 nodupkey;
by dmel6_jxnhash dsim2_jxnHash;
run;


proc sort data=dsim_2_dmel6;
by dmel6_jxnHash;
run;

proc transpose data=dsim_2_dmel6 out=dmel6_with_dsim2 (drop=_name_) prefix=dsim_;
var dsim2_jxnHash;
by dmel6_jxnHash;
run;


data sexsplic.fiveSpecies_dmel6_with_dsim;
merge fivespecies_2dmel_id2 (in=in1) dmel6_with_dsim2(in=in2);
by dmel6_jxnHash;
if in1 then output sexsplic.fiveSpecies_dmel6_with_dsim; 
rename geneid=dmel6_geneid;
keep geneid dmel6_jxnHash flag_dmel650_2_dmel6_ujc flag_dsim202_2_dsim2_ujc flag_dsimwxd_2_dsim2_ujc flag_all mel_symbol flag_sexdet_gene dsim:;
run;


proc export data = sexsplic.fiveSpecies_dmel6_with_dsim
    outfile = "!MCLAB/sex_specific_splicing/cross_species_link_files/fiveSpecies_dmel6_with_dsim.csv" 
    dbms = csv replace ;
run;






