libname sexsplic "!MCLAB/sex_specific_splicing/sasdata";



/*  import fiveSpecies on dsan coord --> has link to dmel jxnHash 
    import link files for all species on the dsan1 genome and the dmel6 genome 
          -- example for dsan11_2_dsan1_ujc_2_dmel6_noGeneID_ujc_xscript_link.csv:  
                has variables dsan1_jxnHash and dmel6_jxnHash

    merge 5species on san to san_2_dmel link file on dsan1_jxnHash
  
    for dsan1_jxnHahess with NO link to dmel6_jxnHash:    
        see if can link via 1) yak, 2) sim2, 3) simw and finally 4) ser
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


    set all together, including remaining dsan1_jxnHashes with no link to a dmel6_jxnhash
    
    create flag for all dsan1_jxnHashs with link to dmel6_jxnHash

    total of 77,804 dsan1_jxnHashes, 12,836 (%16.5) have not link to a dmel6_jxnhash

    *** compared to LMM file - identical!!!  

*/

proc import datafile ="!MCLAB/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_dsan1_ujc_w_IR_flag.csv" 
out =  fivespecies_on_dsan
dbms = csv replace ;
guessingrows = MAX ;
run;

proc sort data = fivespecies_on_dsan ;
by dsan1_jxnHash ;
run;

/* import files linking each species to dsan genome
 
%prep_4_merge (dsan11_2_dsan1_ujc,  dmel6,   dsan1 );
%prep_4_merge (annotation,          genome,  anno  );
 */

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

%import_links (dmel650_2_dmel6_ujc, dsan1,   dmel6 );
%import_links (dser11_2_dser1_ujc,  dsan1,   dser1 );
%import_links (dsim202_2_dsim2_ujc, dsan1,   dsim2 );
%import_links (dsimWXD_2_dsim2_ujc, dsan1,   dsimW );
%import_links (dyak21_2_dyak2_ujc,  dsan1,   dyak2 );


/* merge 5species on san to san_2_dmel link file */

data san_2_mel_linked others_2_track ;
merge fivespecies_on_dsan (in=in1) dsan1_2_dmel6_id (in=in2);
by dsan1_jxnHash;
if in1 and in2 then output san_2_mel_linked;
else  output others_2_track;
run;
    /* 25,220 in san_2_mel_linked
       52,584 in others_2_track */

    
/* others to link:  yak, sim2, simw, mel, and ser */
data yak_onsan_2_mel;
merge dyak2_2_dmel6_id (in=in1) dyak2_2_dsan1_id (in=in2);
by dyak2_jxnHash;
if in1 and in2;
run;  /* 26782 obs */

data sim2_onsan_2_mel;
merge dsim2_2_dmel6_id (in=in1) dsim2_2_dsan1_id (in=in2);
by dsim2_jxnHash;
if in1 and in2;
run; /* 25175 obs */

data simW_onsan_2_mel;
merge dsimW_2_dmel6_id (in=in1) dsimW_2_dsan1_id (in=in2);
by dsimW_jxnHash;
if in1 and in2;
run; /* 25175 obs */

data ser_onsan_2_mel;
merge dser1_2_dmel6_id (in=in1) dser1_2_dsan1_id (in=in2);
by dser1_jxnHash;
if in1 and in2;
run; /* 15508 obs */

/* find where dsan1_jxnHash is not unique ---> keeping first only?????   */
%macro dups (shrt, gen) ;

proc sort data=&shrt._onsan_2_mel ;
by dsan1_jxnHash dmel6_jxnHash &gen._jxnHash;
run;

data &shrt._onsan_2_mel_02 ;
set &shrt._onsan_2_mel ;
by dsan1_jxnHash dmel6_jxnHash &gen._jxnHash;
retain dup_flag 0;
if first.dsan1_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data &shrt._onsan_2_mel_nodup;
set &shrt._onsan_2_mel_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = &shrt._onsan_2_mel_nodup; 
by dsan1_jxnHash ;
run;

%mend ;

%dups (yak, dyak2) ;   /* 26,641 obs */
%dups (sim2, dsim2) ;  /* 24,636 obs */
%dups (simW, dsimW) ;  /* 24,636 obs */
%dups (ser, dser1) ;   /* 14,192 obs */
%dups (mel, dmel6) ;   /* 14,192 obs */

proc sort data=dmel6_2_dsan1_id ;
by dsan1_jxnHash dmel6_jxnHash ;
run;

data dmel6_2_dsan1_id_02 ;
set dmel6_2_dsan1_id ;
by dsan1_jxnHash dmel6_jxnHash ;
retain dup_flag 0;
if first.dsan1_jxnHash then do ;
dup_flag = 0 ;
end;
else do;
dup_flag = 1 ;
end ;
run;

data mel_onsan_2_mel_nodup;
set dmel6_2_dsan1_id_02;
where dup_flag = 0 ;
drop dup_flag ;
run;  

proc sort data = mel_onsan_2_mel_nodup; 
by dsan1_jxnHash ;
run;



/* see if yak, sim2, simW, mel or ser in others_2_track */

proc sort data = others_2_track;
by dsan1_jxnHash ;
run ;

data added_yak others_more ;
merge others_2_track (in=in1) yak_onsan_2_mel_nodup (in=in2) ;
by dsan1_jxnHash ;
if in1 and in2 then output added_yak ;
else if in1 then output others_more ;
run;
    /*  added_yak    10933 obs 
        others_more  41148 obs */

data added_mel others_more_mel ;
merge others_more (in=in1) mel_onsan_2_mel_nodup (in=in2) ;
by dsan1_jxnHash ;
if in1 and in2 then output added_mel ;
else if in1 then output others_more_mel ;
run;
    /*  added_mel        16557 obs 
        others_more_mel  24591 obs */

data added_sim2 others_more3 ;
merge others_more_mel (in=in1) sim2_onsan_2_mel_nodup (in=in2) ;
by dsan1_jxnHash ;
if in1 and in2 then output added_sim2 ;
else if in1 then output others_more3 ;
run;
    /*  added_sim2      9530 obs 
        others_more3   15061 obs */

data added_simW others_more4 ;
merge others_more3 (in=in1) simW_onsan_2_mel_nodup (in=in2) ;
by dsan1_jxnHash ;
if in1 and in2 then output added_simw ;
else if in1 then output others_more4 ;
run;
    /*  added_sim       4778 obs 
        others_more4   10283 obs */

data added_ser nomatch ;
merge others_more4 (in=in1) ser_onsan_2_mel_nodup (in=in2) ;
by dsan1_jxnHash ;
if in1 and in2 then output added_ser ;
else if in1 then output nomatch ;
run;
    /*  added_ser      8485 obs 
        nomatch       1798 obs */

data setting ;
retain dsan1_jxnHash dmel6_jxnHash flag_nomatch ; 
set san_2_mel_linked added_mel added_yak added_sim2 added_simW added_ser nomatch (in=in5) ;
if in5 then flag_nomatch = 1 ;
else flag_nomatch =  0 ;
keep dsan1_jxnHash dmel6_jxnHash flag_nomatch ;
run;

data set2 ;
set setting ;
if flag_nomatch=1 then flag_mel_jxnHash=0;
else flag_mel_jxnHash=1;
drop flag_nomatch ;
run ;

proc sort data=set2 out=test nodupkey;
by dsan1_jxnHash ;
run; /* 77804 obs, no dups */

proc freq data = set2 ;
tables flag_mel_jxnHash ;
run;

data sexsplic.cross_species_dsan_2_dmel6; 
set set2 ;
run ;

title "flag_mel_jxnHash in cross_species_dsan_2_dmel6"; 
proc freq data = sexsplic.cross_species_dsan_2_dmel6; 
tables flag_mel_jxnHash ;
run ;
title "":
/*  flag_mel_jxnHash    Frequency     Percent     Frequency
  --------------------------------------------------------
                 0        1798        2.33          1798
                 1       75225       97.67         77023
*/

/* compare to lmm dataset */
/*
proc freq data = sexsplic.cross_species_dsan_2_dmel6;
tables flag_mel_jxnHash ;
run;


proc sort data = setting ;
by dsan1_jxnHash dmel6_jxnHash flag_mel_jxnHash;
proc sort data = sexsplic.cross_species_dsan_2_dmel6;
by dsan1_jxnHash dmel6_jxnHash flag_mel_jxnHash; 
run ;

proc compare base = test compare = sexsplic.cross_species_dsan_2_dmel6
out = diff_out outnoequal outbase outcompare outdiff ;
run;    /* 0 obs in diff_out!!! */

