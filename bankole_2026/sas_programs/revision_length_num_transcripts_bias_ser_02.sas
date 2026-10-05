

/*from the ser_gene_list_compare code*/


PROC IMPORT OUT= ser_gene
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/gene_summary_dser.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data ser_gene2;
set ser_gene;
rename geneid=dser_geneid;
run;

proc sort data=ser_gene2;
by dser_geneid;
run;

PROC IMPORT OUT= ser_as
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/gene_summary_from_as_analysis_dser.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data ser_as2;
set ser_as;
rename geneid=dser_geneid;
run;


proc sort data=ser_as2;
by dser_geneid;
run;

PROC IMPORT OUT= ser_gene_model
            DATAFILE= "!MCLAB/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dser1_anno_files/dser1_gene_model_info.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;



data ser_gene_model2;
set ser_gene_model;
dser_geneid=geneid;
if num_exonregion le 3 then gene_size="short";
else if num_exonregion le 9 then gene_size="med";
else if num_exonregion >9 then gene_size="long";

proc sort data=ser_gene_model2;
by dser_geneid;
run;


data ser_results;
merge ser_gene2 (in=in1) ser_gene_model2  ser_as2;
by dser_geneid;
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

proc freq data=ser_results;
where sexbias ne "not_evaluated" and lots ge 1;
tables lots*gene_size/chisq cmh;
run;
/*p<0.0001 long genes those with  more exons have more transcripts;*/

proc freq data=ser_results;
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

proc freq data=ser_results;
where lots=3;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;
/*this is significant p=0.4557, 0.2506,0.1473,.2824,0.7140,0.0458*/

proc freq data=ser_results;
where lots=2;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;

/*not significant alt 3p 0.0236*/


proc freq data=ser_results;
where lots=1;
tables sexBias1*(flag_alt_dataonly 
flag_alt_ir
flag_alt_donor_acceptor
flag_alt_5p_er
flag_alt_3p_er
flag_alt_erskip )/chisq;
run;

/* p=0.0001, 0.0001, 0.0001, 0.0560,0.0012,0.0001*/

proc freq data=ser_results;
*where sexbias ne "not_evaluated" and lots=1;
*where sexbias ne "not_evaluated" and lots=2;
*where sexbias ne "not_evaluated" and lots=3;
*where sexbias ne "not_evaluated" and lots>1;
tables sexbias1*gm_size/chisq;
run;

/*lots=1,2,3 p=0.0011, 0.0098 0.0001, >1  0.0001*/

                  
proc freq data=ser_results;
where sexbias ne "not_evaluated" and lots>0;

tables lots*sexBias1*gene_size/cmh chisq;
run;                        

proc freq data=ser_results;
where lots >2;
tables sexbias1*gene_size/chisq;
run;

/*p<0.0001*/

     
                   
proc freq data=ser_results;
where sexbias ne "not_evaluated" and lots>0 and sexbias1 ne "O";

tables lots*sexBias1*gene_size/cmh chisq;
run;            



/*alt data only signficant*/

sexBias1     lots

                                         Frequency|
                                         Percent  |
                                         Row Pct  |
                                         Col Pct  |       1|       2|       3|  Total
                                         ---------+--------+--------+--------+
                                         B        |     66 |    144 |    397 |    607
                                                  |   0.75 |   1.63 |   4.50 |   6.89
                                                  |  10.87 |  23.72 |  65.40 |
                                                  |   1.51 |   5.77 |  20.23 |
                                         ---------+--------+--------+--------+
                                         O        |   1454 |   1353 |   1265 |   4072
                                                  |  16.49 |  15.35 |  14.35 |  46.19
                                                  |  35.71 |  33.23 |  31.07 |
                                                  |  33.36 |  54.21 |  64.48 |
                                         ---------+--------+--------+--------+
                                         unbiased |   2838 |    999 |    300 |   4137
                                                  |  32.19 |  11.33 |   3.40 |  46.93
                                                  |  68.60 |  24.15 |   7.25 |
                                                  |  65.12 |  40.02 |  15.29 |
                                         ---------+--------+--------+--------+
                                         Total        4358     2496     1962     8816
                                                     49.43    28.31    22.25   100.00




proc freq data=ser_results;
where sexbias ne "not_evaluated";
tables sexbias*gene_size;
tables sexBias1*gene_size/cmh chisq;
run;
 				/*	sexBias        gene_size more than 9, 4-9, 3 or fewer */ 

                                     sexBias1     gene_size

                                         Frequency|
                                         Percent  |
                                         Row Pct  |
                                         Col Pct  |long    |med     |short   |  Total
                                         ---------+--------+--------+--------+
                                         B        |    220 |    314 |     74 |    608
                                                  |   1.93 |   2.76 |   0.65 |   5.35
                                                  |  36.18 |  51.64 |  12.17 |
                                                  |  13.77 |   6.38 |   1.52 |
                                         ---------+--------+--------+--------+
                                         O        |    810 |   2171 |   1429 |   4410
                                                  |   7.12 |  19.09 |  12.57 |  38.78
                                                  |  18.37 |  49.23 |  32.40 |
                                                  |  50.69 |  44.14 |  29.43 |
                                         ---------+--------+--------+--------+
                                         unbiased |    568 |   2433 |   3352 |   6353
                                                  |   5.00 |  21.40 |  29.48 |  55.87
                                                  |   8.94 |  38.30 |  52.76 |
                                                  |  35.54 |  49.47 |  69.04 |
                                         ---------+--------+--------+--------+
                                         Total        1598     4918     4855    11371
                                                     14.05    43.25    42.70   100.00

proc freq data=ser_results;
where sexbias ne "not_evaluated" and lots=3;

tables sexBias1*gene_size/cmh chisq;
run;                              
                    
               *     p<0.0001; lots >1 p<0.0001, lots =3 p=0.0001;
                    
                    
proc freq data=ser_results;                
where sexbias ne "not_evaluated" and lots>0;

tables lots*sexBias2*gene_size/cmh chisq;
run;      
/*row mean scores 0.001,  for strata 1,3 not 2*/

proc sort data=ser_results;
by gm_size;

proc univariate data=ser_results normal plot;  
by gm_size;
var num_ujc_analyzable;
run;
                        
                        
                        proc plot data=ser_results;
                        where num_ujc_analyzable <20 and num_exonregion <20;
                        plot num_ujc_analyzable*num_exonregion;
                        run;
                                 
                                                           
proc sort data=ser_results;
by sexbias;

proc freq data=ser_results;
where num_exonregion >9;
tables sexbias1*(alt_dataonly alt_ir alt_donor_acceptor alt_5p_er alt_3p_er alt_erskip)/chisq;
run;

proc univariate data=ser_results normal plot;
where sexbias ne "not_evaluated";
by sexbias;
var num_exonregion;
run;


proc univariate data=ser_results normal plot;
where sexbias ne "not_evaluated";
by sexbias;
var num_unique_erp;
run;




