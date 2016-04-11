/* Checking GATK output, compare LD and TD
Need to prepare the files a little before plotting */

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

*Import TD and LD files;

 data WORK.TD_x    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/checking_gatk_output/tajima_d_gvcf.txt' delimiter='09'x
  MISSOVER DSD lrecl=32767 ;
         informat chrom $1. ;
         informat bin_start best32. ;
         informat n_snps best32. ;
         informat TajimaD best32. ;
         format chrom $1. ;
         format bin_start best12. ;
         format n_snps best12. ;
         format TajimaD best12. ;
      input
                  chrom $
                  bin_start
                  n_snps
                  TajimaD
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

 data WORK.LD_x    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/checking_gatk_output/chromX_ld_window_avg.txt' delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat chrom $1. ;
        informat bin best32. ;
        informat R2 best32. ;
        format chrom $1. ;
        format bin best12. ;
        format R2 best12. ;
     input
                 chrom $
                 bin $
                 R2
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data td1;
  set td_x;
  bin = bin_start/10000;
  run;

proc sort data=ld_x;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td_x;
  merge td1 (in=in1) ld_x (in=in2);
  by bin;
  run;

proc export data=ld_td_x
    outfile="/home/fnew/dsim/checking_gatk_output/chrX_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;

