/* Manipulate the tajima's d file. First for 100kb window*/

libname dsim "/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/sas_data";


data WORK.tsd100    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile ‘!MCLAB/ethanol/Sim_Pop_Gen/output/tajima_d/filter_nolab_tajd_100kb.Tajima.D' delimiter='09'x
 MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat CHROM $2. ;
         informat BIN_START best32. ;
         informat N_SNPS best32. ;
         informat TajimaD best32. ;
         format CHROM $2. ;
         format BIN_START best12. ;
         format N_SNPS best12. ;
         format TajimaD best12. ;
      input

                  CHROM $
                  BIN_START
                  N_SNPS
                  TajimaD
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


data WORK.tsd10    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '!MCLAB/ethanol/Sim_Pop_Gen/output/tajima_d/filter_nolab_tajd_10kb.Tajima.D' delimiter='09'x
 MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat CHROM $2. ;
         informat BIN_START best32. ;
         informat N_SNPS best32. ;
         informat TajimaD best32. ;
         format CHROM $2. ;
         format BIN_START best12. ;
         format N_SNPS best12. ;
         format TajimaD best12. ;
      input

                  CHROM $
                  BIN_START
                  N_SNPS
                  TajimaD
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


*Set significance levels. flag if sig;
*n=175, critical values: -1.765~2.095;

data tsd100_sig;
    set tsd100;
    if TajimaD < -1.765 then flag_sig=1; else if TajimaD > 2.095 then flag_sig=1;
        else flag_sig=0;
    run;

proc freq data=tsd100_sig;
  tables flag_sig;
  run; *1165 nonsig, 22 sig;


data tsd10_sig;
    set tsd10;
    if TajimaD < -1.765 then flag_sig=1; else if TajimaD > 2.095 then flag_sig=1;
        else flag_sig=0;
    run;

proc freq data=tsd10_sig;
  tables flag_sig;
  run; *10767 nonsig, 1080 sig;

  proc freq data=tsd10_sig;
    tables CHROM*flag_sig;
    run;

proc sort data=tsd10_sig;
  by TajimaD;
  run; *range: -2.906 to 4.44;

proc sort data=tsd100_sig;
  by TajimaD;
  run; *range: -2.903 to 3.12;

data tsd100_bin;
  set tsd100_sig;
  if TajimaD > -0.1 and TajimaD < 0.1 then bin = 0;
  else if TajimaD > -1 and TajimaD <= -0.1 then bin = -1;
  else if TajimaD >-1.765 and TajimaD <= -1 then bin= -2;
  else if TajimaD > -3 and TajimaD <= -1.765 then bin = -3; *sig;
  else if TajimaD >= 0.1 and TajimaD < 1 then bin = 1;
  else if TajimaD >=1 and TajimaD < 2.095 then bin=2;
  else if TajimaD >=2.095 and TajimaD < 3 then bin=3; *sig;
  else if TajimaD >=3 and TajimaD < 4 then bin=4; 
  else if TajimaD >=4 and TajimaD < 5 then bin=5;
  run;

data tsd10_bin;
  set tsd10_sig;
  if TajimaD > -0.1 and TajimaD < 0.1 then bin = 0;
  else if TajimaD > -1 and TajimaD <= -0.1 then bin = -1;
  else if TajimaD >-1.765 and TajimaD <= -1 then bin= -2;
  else if TajimaD > -3 and TajimaD <= -1.765 then bin = -3; *sig;
  else if TajimaD >= 0.1 and TajimaD < 1 then bin = 1;
  else if TajimaD >=1 and TajimaD < 2.095 then bin=2;
  else if TajimaD >=2.095 and TajimaD < 3 then bin=3; *sig;
  else if TajimaD >=3 and TajimaD < 4 then bin=4; 
  else if TajimaD >=4 and TajimaD < 5 then bin=5;
  run;

proc freq data = tsd100_bin;
    tables flag_sig*CHROM;
    run;
proc freq data=tsd10_bin;
    tables flag_sig*CHROM;
    run;

proc freq data = tsd100_bin;
    tables bin*CHROM;
    run;
proc freq data=tsd10_bin;
    tables bin*CHROM;
    run;
