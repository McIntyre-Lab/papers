

/*from the mel_gene_list_compare code*/



PROC IMPORT OUT= mel_as
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/gene_summary_from_as_analysis_dmel.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data mel_as2;
set mel_as;
rename geneid=dmel_geneid;
run;


PROC IMPORT OUT= mel_gene_model
            DATAFILE= "!MCLAB/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dmel6_anno_files/dmel6_gene_model_info.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

proc contents data=mel_gene_model;
run;

proc univariate data=mel_gene_model normal plot ;
var num_exonregion;
run;


data mel_gene_model2;
set mel_gene_model;
dmel_geneid=geneid;
if num_exonregion le 3 then gene_size="short";
else if num_exonregion le 9 then gene_size="med";
else if num_exonregion >9 then gene_size="long";


proc univariate data=mel_gene2 normal plot ;

var num_ujc_analyzable;
run;


proc sort data=mel_gene_model2;
by dmel_geneid;

data mel_results;
merge mel_gene2 (in=in1) mel_gene_model2  mel_as2;
by dmel_geneid;
if in1;

flag_alt_dataonly =.;
flag_alt_ir =.;
flag_alt_donor_acceptor =.;
flag_alt_5p_er=.;
flag_alt_3p_er =.;
flag_alt_erskip=.;

if sexBias="unbiased" then sexBias1="unbiased";
if sexBias="B" then sexBias1="B";
if sexBias="M" or sexBias="F" then sexBias1="O";


if sexBias="unbiased" then sexBias2="unbiased";
if sexBias="B" then sexBias2="B";
if sexBias="M" or sexBias="F" then sexBias2="B";


if numExon_gm le 3 then gm_size="short";
else if numExon_gm le 9 then gm_size="med";
else if numExon_gm  >9 then gm_size="long";

if num_ujc_analyzable >8 then lots=3;
else if num_ujc_analyzable >4 then lots=2;
else if num_ujc_analyzable >1 then lots=1;
else if num_ujc_analyzable le 1 then lots=0;


if alt_dataonly ="none" then flag_alt_dataonly  =0; else if alt_dataonly  ne "." then flag_alt_dataonly =1;
if alt_ir="none"  then flag_alt_ir =0; else if alt_ir ne "." then flag_alt_ir=1;
if alt_donor_acceptor ="none" then flag_alt_donor_acceptor  =0; else if alt_donor_acceptor  ne "." then flag_alt_donor_acceptor =1;
if alt_5p_er ="none" then flag_alt_5p_er=0; else if alt_5p_er ne "." then flag_alt_5p_er =1;
if alt_3p_er ="none" then flag_alt_3p_er =0; else if alt_3p_er ne "." then flag_alt_3p_er=1;
if alt_erskip ="none"  then flag_alt_erskip  =0; else if alt_erskip  ne "." then flag_alt_erskip =1;


run;

proc freq data=mel_results;
where sexbias ne "not_evaluated" and lots ge 1;
tables lots*gene_size/chisq cmh;
run;
/*p<0.0001 long genes have more exons;*/

proc freq data=mel_results;
where sexbias ne "not_evaluated" and lots ge 1;
tables sexBias1*lots;
tables lots
*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;

proc freq data=mel_results;
where lots=3;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;
/*nothing is significnat here*/

proc freq data=mel_results;
where lots=2;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;

/*alt 5 p=0.0026 and alt3 p=0.0071*/


proc freq data=mel_results;
where lots=1;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;
/*0.05, 0.1453, 0.2447, 0.1464,0.2082,0.002*/

proc freq data=mel_results;
where lots=2;
tables sexbias1*gm_size/chisq;
run;

/*p=0.0003*/

/*alt data only signficant*/

proc freq data=mel_results;
where lots >0;
tables sexbias1*lots/chisq;
run;

/*p<0.0001*/


 sexBias1     lots

                                         Frequency|
                                         Percent  |
                                         Row Pct  |
                                         Col Pct  |       1|       2|       3|  Total
                                         ---------+--------+--------+--------+
                                         B        |     44 |    107 |    257 |    408
                                                  |   0.51 |   1.24 |   2.98 |   4.73
                                                  |  10.78 |  26.23 |  62.99 |
                                                  |   0.98 |   4.29 |  15.73 |
                                         ---------+--------+--------+--------+
                                         O        |   1428 |   1247 |   1081 |   3756
                                                  |  16.54 |  14.44 |  12.52 |  43.51
                                                  |  38.02 |  33.20 |  28.78 |
                                                  |  31.69 |  50.02 |  66.16 |
                                         ---------+--------+--------+--------+
                                         unbiased |   3034 |   1139 |    296 |   4469
                                                  |  35.14 |  13.19 |   3.43 |  51.77
                                                  |  67.89 |  25.49 |   6.62 |
                                                  |  67.33 |  45.69 |  18.12 |
                                         ---------+--------+--------+--------+
                                         Total        4506     2493     1634     8633
                                                     52.20    28.88    18.93   100.00



proc freq data=mel_results;
tables gm_size*gene_size;
run;



proc freq data=mel_results;
where sexbias ne "not_evaluated";
tables sexbias*gene_size;
tables sexBias1*gene_size/cmh chisq;
run;
 				/*	sexBias        gene_size more than 7, 4-7, 3 or fewer */ 

                                      Frequency     |
                                      Percent       |
                                      Row Pct       |
                                      Col Pct       |long    |med     |short   |  Total
                                      --------------+--------+--------+--------+
                                      B             |    204 |    143 |     63 |    410
                                                    |   1.80 |   1.26 |   0.56 |   3.62
                                                    |  49.76 |  34.88 |  15.37 |
                                                    |   8.29 |   3.56 |   1.30 |
                                      --------------+--------+--------+--------+
                                      F             |    361 |    837 |    829 |   2027
                                                    |   3.19 |   7.40 |   7.33 |  17.92
                                                    |  17.81 |  41.29 |  40.90 |
                                                    |  14.66 |  20.86 |  17.14 |
                                      --------------+--------+--------+--------+
                                      M             |    761 |    746 |    516 |   2023
                                                    |   6.73 |   6.59 |   4.56 |  17.88
                                                    |  37.62 |  36.88 |  25.51 |
                                                    |  30.91 |  18.59 |  10.67 |
                                      --------------+--------+--------+--------+
                                      unbiased      |   1136 |   2287 |   3429 |   6852
                                                    |  10.04 |  20.22 |  30.31 |  60.57
                                                    |  16.58 |  33.38 |  50.04 |
                                                    |  46.14 |  56.99 |  70.89 |
                                      --------------+--------+--------+--------+
                                      Total             2462     4013     4837    11312
                                                       21.76    35.48    42.76   100.00          sexBias1     gene_size

/*now with 9 in the deifnition*/
                                         Frequency|
                                         Percent  |
                                         Row Pct  |
                                         Col Pct  |long    |med     |short   |  Total
                                         ---------+--------+--------+--------+
                                         B        |    143 |    204 |     63 |    410
                                                  |   1.26 |   1.80 |   0.56 |   3.62
                                                  |  34.88 |  49.76 |  15.37 |
                                                  |   9.15 |   4.15 |   1.30 |
                                         ---------+--------+--------+--------+
                                         O        |    741 |   1964 |   1345 |   4050
                                                  |   6.55 |  17.36 |  11.89 |  35.80
                                                  |  18.30 |  48.49 |  33.21 |
                                                  |  47.44 |  39.98 |  27.81 |
                                         ---------+--------+--------+--------+
                                         unbiased |    678 |   2745 |   3429 |   6852
                                                  |   5.99 |  24.27 |  30.31 |  60.57
                                                  |   9.89 |  40.06 |  50.04 |
                                                  |  43.41 |  55.87 |  70.89 |
                                         ---------+--------+--------+--------+
                                         Total        1562     4913     4837    11312
                                                     13.81    43.43    42.76   100.00
                                                       
                         
proc freq data=mel_results;
where sexbias ne "not_evaluated" and lots>1;

tables sexBias1*gene_size/cmh chisq;
run;                              
                    
               *     p<0.0001; lots >1 p<0.0001, lots =3 p=0.06;
                    
    
proc freq data=mel_results;                
where sexbias ne "not_evaluated" and lots>0;

tables lots*sexBias1*gene_size/cmh chisq;
run;                        
                    /*row mean scores p=.0023  significant when contorlling for lots=2 but not for 1,3*/
                
                
proc freq data=mel_results;
where lots >3;
tables sexbias1*gene_size/chisq;
run;

/*p<0.0001*/


proc freq data=mel_results;
where sexbias ne "not_evaluated" and lots>0 and sexbias1 ne "O";

tables lots*sexBias1*gene_size/cmh chisq;
run;            

/*not signficnt*/

proc freq data=mel_results;
where sexbias ne "not_evaluated" and lots>0 ;

tables lots*sexBias2*gene_size/cmh chisq;
run;            

  
proc sort data=mel_results;
by gm_size;

proc univariate data=mel_results normal plot;  
by gm_size;
var num_ujc_analyzable;
run;
                        
                        
                        proc plot data=mel_results;
                        where num_ujc_analyzable <20 and num_exonregion <20;
                        plot num_ujc_analyzable*num_exonregion;
                        run;
                                 
                                                           
proc sort data=mel_results;
by sexbias;

proc freq data=mel_results;
where num_exonregion >9;
tables sexbias1*(alt_dataonly alt_ir alt_donor_acceptor alt_5p_er alt_3p_er alt_erskip)/chisq;
run;

proc univariate data=mel_results normal plot;
where sexbias ne "not_evaluated";
by sexbias;
var num_exonregion;
run;


proc univariate data=mel_results normal plot;
where sexbias ne "not_evaluated";
by sexbias;
var num_unique_erp;
run;




