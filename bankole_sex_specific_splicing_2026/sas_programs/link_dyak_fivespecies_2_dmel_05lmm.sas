libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";



/*  import fiveSpecies on dyak coord --> has link to dmel jxnHash 
    import link files for all species on the dyak2 genome and the dmel6 genome 
          -- example for dsan11_2_dsan1_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv:  
                has variables dsan1_jxnHash and dmel6_jxnHash

    merge 5species on yak to yak_2_dmel link file on dyak_jxnHash
  
    for dyak_jxnHahess with NO link to dmel6_jxnHash:    
        see if can link via 1) san, 2) sim2, 3) simw and finally 4) ser
        *** note that order of merges is important! by descending relatedness!

               +----------- D. yakuba
               |
      +--------+           
      |        +----------- D. santomea
      |
      |                +--- D. simulans
      |        +-------+
      |        |       +--- D. melanogaster
      +--------+
               |
               +----------- D. serrata


    set all together, including remaining dyak_jxnHashes with no link to a dmel6_jxnhash
    
    create flag for all dyak_jxnHashs with link to dmel6_jxnHash

                                                 Cumulative
    flag_mel_jxnHash    Frequency     Percent     Frequency
    --------------------------------------------------------
                   0       12807       16.45         12807
                   1       65052       83.55         77859
   


*/

proc import datafile ="!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_dyak2_ujc_w_IR_flag.csv" 
out =  fivespecies_on_dyak
dbms = csv replace ;
guessingrows = MAX ;
run;

proc sort data = fivespecies_on_dyak ;
by dyak2_jxnHash ;
run;

/* import files linking each species to dyak genome and dmel6 genome*/

%macro import_links (annotation,          genome,  anno  );

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
run;

proc sort data=&anno._2_&genome._id;
by &anno._jxnHash;
run;

%mend ;

%import_links (dsan11_2_dsan1_ujc,  dmel6,   dsan1 );
%import_links (dser11_2_dser1_ujc,  dmel6,   dser1 );
%import_links (dsim202_2_dsim2_ujc, dmel6,   dsim2 );
%import_links (dsimWXD_2_dsim2_ujc, dmel6,   dsimW );
%import_links (dyak21_2_dyak2_ujc,  dmel6,   dyak2 );

%import_links (dmel650_2_dmel6_ujc, dyak2,   dmel6 );
%import_links (dser11_2_dser1_ujc,  dyak2,   dser1 );
%import_links (dsim202_2_dsim2_ujc, dyak2,   dsim2 );
%import_links (dsimWXD_2_dsim2_ujc, dyak2,   dsimW );
%import_links (dsan11_2_dsan1_ujc,  dyak2,   dsan1 );


/* merge 5species on san to san_2_dmel link file */

data yak_2_mel_linked others_2_track ;
merge fivespecies_on_dyak (in=in1) dyak2_2_dmel6_id (in=in2);
by dyak2_jxnHash;
if in1 and in2 then output yak_2_mel_linked;
else  output others_2_track;
run;
    /* 26,908 in yak_2_mel_linked
       50,951 in others_2_track */

    
/* others to link: 1) san, 2) sim2, 3) simw and finally 4) ser  */
data san_onyak_2_mel;
merge dsan1_2_dmel6_id (in=in1) dsan1_2_dyak2_id (in=in2);
by dsan1_jxnHash;
if in1 and in2;
run; /* 25214 obs */

data sim2_onyak_2_mel;
merge dsim2_2_dmel6_id (in=in1) dsim2_2_dyak2_id (in=in2);
by dsim2_jxnHash;
if in1 and in2;
run; /* 25190 obs */

data simW_onyak_2_mel;
merge dsimW_2_dmel6_id (in=in1) dsimW_2_dyak2_id (in=in2);
by dsimW_jxnHash;
if in1 and in2;
run; /* 21273 obs */

data ser_onyak_2_mel;
merge dser1_2_dmel6_id (in=in1) dser1_2_dyak2_id (in=in2);
by dser1_jxnHash;
if in1 and in2;
run;  /* 15514 obs */



/* find where dsan1_jxnHash is not unique ---> keeping first only?????   */
%macro dups (shrt, gen) ;

proc sort data=&shrt._onyak_2_mel ;
by dyak2_jxnHash dmel6_jxnHash &gen._jxnHash;
run;

data &shrt._onyak_2_mel_02 ;
set &shrt._onyak_2_mel ;
by dyak2_jxnHash dmel6_jxnHash &gen._jxnHash;
retain dup_flag 0;
if first.dyak2_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data &shrt._onyak_2_mel_nodup;
set &shrt._onyak_2_mel_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = &shrt._onyak_2_mel_nodup; 
by dyak2_jxnHash ;
run;

%mend ;

%dups (san, dsan1) ;   /* 25,001 obs */
%dups (sim2, dsim2) ;  /* 24,667 obs */
%dups (simW, dsimW) ;  /* 21059 obs */
%dups (ser, dser1) ;   /* 14,211 obs */


proc sort data=dmel6_2_dyak2_id ;
by dyak2_jxnHash dmel6_jxnHash ;
run;

data dmel6_2_dyak2_id_02 ;
set dmel6_2_dyak2_id ;
by dyak2_jxnHash dmel6_jxnHash ;
retain dup_flag 0;
if first.dyak2_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data mel_onyak_2_mel_nodup;
set dmel6_2_dyak2_id_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = mel_onyak_2_mel_nodup; 
by dyak2_jxnHash ;
run;



/* see if 1) san, 2) sim2, 3) simw and finally 4) ser in others_2_track */

proc sort data = others_2_track;
by dyak2_jxnHash ;
run ;

data added_san others_more ;
merge others_2_track (in=in1) san_onyak_2_mel_nodup (in=in2) ;
by dyak2_jxnHash ;
if in1 and in2 then output added_san ;
else if in1 then output others_more ;
run;
    /*  added_san    9351 obs 
        others_more  41096 obs */

data added_mel others_more_mel ;
merge others_more (in=in1) mel_onyak_2_mel_nodup (in=in2) ;
by dyak2_jxnHash ;
if in1 and in2 then output added_mel ;
else if in1 then output others_more_mel ;
run;
    /*  added_mel        16554 obs 
        others_more_mel  24542 obs */

data added_sim2 others_more3 ;
merge others_more_mel (in=in1) sim2_onyak_2_mel_nodup (in=in2) ;
by dyak2_jxnHash ;
if in1 and in2 then output added_sim2 ;
else if in1 then output others_more3 ;
run;
    /*  added_sim2      9482 obs 
        others_more3    15060 obs */

data added_simW others_more4 ;
merge others_more3 (in=in1) simW_onyak_2_mel_nodup (in=in2) ;
by dyak2_jxnHash ;
if in1 and in2 then output added_simw ;
else if in1 then output others_more4 ;
run;
    /*  added_simw     4756 obs 
        others_more4   10304 obs */

data added_ser nomatch ;
merge others_more4 (in=in1) ser_onyak_2_mel_nodup (in=in2) ;
by dyak2_jxnHash ;
if in1 and in2 then output added_ser ;
else if in1 then output nomatch ;
run;
    /*  added_ser     8523 obs 
        nomatch       1781 obs */

data setting ;
retain dyak2_jxnHash dmel6_jxnHash flag_nomatch ; 
set yak_2_mel_linked added_mel added_san added_sim2 added_simW added_ser nomatch (in=in5) ;
if in5 then flag_nomatch = 1 ;
else flag_nomatch =  0 ;
keep dyak2_jxnHash dmel6_jxnHash flag_nomatch ;
run;

data set2 ;
set setting ;
if flag_nomatch=1 then flag_mel_jxnHash=0;
else flag_mel_jxnHash=1;
drop flag_nomatch ;
run ;

proc sort data=set2 out=test nodupkey;
by dyak2_jxnHash ;
run; /* 77804 obs, no dups */

proc freq data = set2 ;
tables flag_mel_jxnHash ;
run;
/*                                  Cumulative
flag_mel_jxnHash    Frequency         Frequency
-------------------------------------------------------
               0       1781               1781
               1       75259              77040
*/

data sexsplic.cross_species_dyak_2_dmel6  ;
set set2 ;
run;

title "flag_mel_jxnHash in cross_species_dyak_2_dmel6"; 
proc freq data = sexsplic.cross_species_dyak_2_dmel6 ;
tables flag_mel_jxnHash ;
run;
title "";

/*added this to export the file for plotting*/

proc export data = sexsplic.cross_species_dyak_2_dmel6 
outfile = "/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/cross_species_dyak_2_dmel6.csv"
dbms = csv replace ;
run;



