libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";



/*  import fiveSpecies on dsim2 coord --> has link to dmel jxnHash 
    import link files for all species on the dsim2 genome and the dmel6 genome 
          -- example for dser1_2_dser1_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv:  
                has variables dser1_jxnHash and dmel6_jxnHash

    merge 5species on sim2 to sim_2_dmel link file on dsim2_jxnHash
  
    for dsim2_jxnHashes with NO link to dmel6_jxnHash:    
        see if can link via 1) ser, 2) san and finally 3) yak  [ no link file for dsimW to dsim2!! ]
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


    set all together, including remaining dsim2_jxnHashes with no link to a dmel6_jxnhash
    
    create flag for all dsim2_jxnHashs with link to dmel6_jxnHash

                                               Cumulative
 flag_mel_jxnHash    Frequency     Percent     Frequency
 -------------------------------------------------------
                0       11558       14.85         11558
                1       66250       85.15         77808                                          
*/

proc import datafile ="!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_dsim2_anno_files/flag_fiveSpecies_2_dsim2_ujc_w_IR_flag.csv" 
out =  fivespecies_on_dsim2
dbms = csv replace ;
guessingrows = MAX ;
run;

proc sort data = fivespecies_on_dsim2 ;
by dsim2_jxnHash ;
run;


proc freq data=fivespecies_on_dsim2;
tables flag_dsimwxd_2_dsim2_ujc flag_dsim202_2_dsim2_ujc flag_dmel650_2_dmel6_ujc flag_dyak21_2_dyak2_ujc flag_dser11_2_dser1_ujc flag_dsan11_2_dsan1_ujc ;
run;



/* import files linking each species to genome */

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

%import_links (dmel650_2_dmel6_ujc, dsim2,   dmel6 );
%import_links (dsan11_2_dsan1_ujc,  dsim2,   dsan1 );
%import_links (dser11_2_dser1_ujc,  dsim2,   dser1 );
%import_links (dyak21_2_dyak2_ujc,  dsim2,   dyak2 );

/* dsimw on dsim2 genome */
data dsimw2_2_dmel6_id ;
set  dsimw_2_dmel6_id ;
rename dsimW_jxnHash = dsim2_jxnHash ;
run;

/* merge 5species on sim2 to the sim_2_dmel link file ==> this should pick up both dsim2 AND dsimW since they are both on dsim2 coords */
data sim_2_mel_linked others_2_track ;
merge fivespecies_on_dsim2 (in=in1) dsim2_2_dmel6_id (in=in2) dsimw2_2_dmel6_id (in=in3) ;
by dsim2_jxnHash;
if (in1 and in2) or (in1 and in3) then output sim_2_mel_linked;
else output others_2_track ;
run;  /* 32322 */


/* prep others to link: 1) ser, 3) san and finally 4) yak */
data ser_onsim_2_mel;
merge dser1_2_dmel6_id (in=in1) dser1_2_dsim2_id (in=in2);
by dser1_jxnHash;
if in1 and in2;
run; /* 15496 obs */

data san_onsim_2_mel;
merge dsan1_2_dmel6_id (in=in1) dsan1_2_dsim2_id (in=in2);
by dsan1_jxnHash;
if in1 and in2;
run; /* 24792 obs */

data yak_onsim_2_mel;
merge dyak2_2_dmel6_id (in=in1) dyak2_2_dsim2_id (in=in2);
by dyak2_jxnHash;
if in1 and in2;
run;  /* 26402 obs */


/* find where dsan1_jxnHash is not unique ---> keeping first only?????   */
%macro dups (shrt, gen) ;

proc sort data=&shrt._onsim_2_mel ;
by dsim2_jxnHash dmel6_jxnHash &gen._jxnHash;
run;

data &shrt._onsim_2_mel_02 ;
set &shrt._onsim_2_mel ;
by dsim2_jxnHash dmel6_jxnHash &gen._jxnHash;
retain dup_flag 0;
if first.dsim2_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data &shrt._onsim_2_mel_nodup;
set &shrt._onsim_2_mel_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = &shrt._onsim_2_mel_nodup; 
by dsim2_jxnHash ;
run;

%mend ;

%dups (yak, dyak2) ;   /* 25,720 obs */
%dups (ser, dser1) ;   /* 14,192 obs */
%dups (san, dsan1) ;   /* 24,200 obs */

proc sort data=dmel6_2_dsim2_id ;
by dsim2_jxnHash dmel6_jxnHash ;
run;

data dmel6_2_dsim2_id_02 ;
set dmel6_2_dsim2_id ;
by dsim2_jxnHash dmel6_jxnHash ;
retain dup_flag 0;
if first.dsim2_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data dmel_onsim_2_mel_nodup;
set dmel6_2_dsim2_id_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = dmel_onsim_2_mel_nodup; 
by dsim2_jxnHash ;
run;


/* see if  1) mel 1.5) ser 2) san 3) yak in others_2_track */

proc sort data = others_2_track;
by dsim2_jxnHash ;
run ;

/* trying adding mel */
data added_mel others_more_w_mel ;
merge others_2_track (in=in1) dmel_onsim_2_mel_nodup (in=in2) ;
by dsim2_jxnHash ;
if in1 and in2 then output added_mel ;
else if in1 then output others_more_w_mel ;
run;
    /*  added_mel           11524 obs 
        others_more_w_mel   30435 obs */


data added_ser others_more ;
merge others_more_w_mel (in=in1) ser_onsim_2_mel_nodup (in=in2) ;
by dsim2_jxnHash ;
if in1 and in2 then output added_ser ;
else if in1 then output others_more ;
run;
    /*  added_ser      9212 obs 
        others_more   21223 obs */

data added_san others_more2 ;
merge others_more (in=in1) san_onsim_2_mel_nodup (in=in2) ;
by dsim2_jxnHash ;
if in1 and in2 then output added_san ;
else if in1 then output others_more2 ;
run;
    /*  added_san      10944 obs 
        others_more2   10279 obs */

data added_yak nomatch ;
merge others_more2 (in=in1) yak_onsim_2_mel_nodup (in=in2) ;
by dsim2_jxnHash ;
if in1 and in2 then output added_yak ;
else if in1 then output nomatch ;
run;
    /*  added_yak    9339 obs 
        nomatch     940 obs */

data setting ;
retain dsim2_jxnHash dmel6_jxnHash flag_nomatch ; 
set sim_2_mel_linked added_mel added_ser added_san added_yak nomatch (in=in5) ;
if in5 then flag_nomatch = 1 ;
else flag_nomatch =  0 ;
keep dsim2_jxnHash dmel6_jxnHash flag_nomatch ;
run;

data set2 ;
set setting ;
if flag_nomatch=1 then flag_mel_jxnHash=0;
else flag_mel_jxnHash=1;
drop flag_nomatch ;
run ;

proc sort data=set2 out=test nodupkey;
by dsim2_jxnHash ;
run; /* 77,281 obs, no dups */

proc freq data = set2 ;
tables flag_mel_jxnHash ;
run; 

data sexsplic.cross_species_dsim2_2_dmel6 ; 
set set2 ;
run ;

/*                                             Cumulative
 flag_mel_jxnHash    Frequency          Frequency
 -------------------------------------------------------
                0       940                940
                1       76341              77281
                                          
*/
title "flag_mel_jxnHash in cross_species_dsim2_2_dmel6"; 
proc freq data = sexsplic.cross_species_dsim2_2_dmel6 ;
tables flag_mel_jxnHash ;
run;
title "":



