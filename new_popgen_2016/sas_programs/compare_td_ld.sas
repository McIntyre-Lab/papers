libname dsim "/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/sas_data";

/* The files “all_chrom#_tajimaD_10kb.Tajima.D are located on the HPC: dsim_pop/vcftools_output/tajima_d_split_vcf/all_chrom*_tajimaD_10kb.Tajima.D */

/* The files “chrom*_ld_window_avg.txt” are located on the HPC: dsim_pop/plink/split_vcf_ld/chrom*_ld_window_avg.txt */


/* chromosome 4*/

   data WORK.LD    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/chrom4_ld_window_avg.txt' delimiter = ',' MISSOVER
 DSD lrecl=32767 ;
        informat size $4. ;
        informat chrom4 $6. ;
        informat bin best32. ;
        informat ld best32. ;
        format size $4. ;
        format chrom4 $6. ;
        format bin best12. ;
        format ld best12. ;
     input
                 size $
                 chrom4 $
                 bin
                 ld
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


  data WORK.TD    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/all_chrom4_10kb_window_tsd.txt' delimiter='09'x
 MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat CHROM best32. ;
        informat BIN_START best32. ;
        informat N_SNPS best32. ;
        informat TajimaD best32. ;
        format CHROM best12. ;
        format BIN_START best12. ;
        format N_SNPS best12. ;
        format TajimaD best12. ;
     input
                 CHROM
                 BIN_START
                 N_SNPS
                 TajimaD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data td1;
  set td;
  bin = bin_start/10000;
  run;

proc sort data=ld;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td;
  merge td1 (in=in1) ld (in=in2);
  by bin;
  run;

proc export data=ld_td
    outfile=“!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chr4_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;


/* Chrom 2L */

   data WORK.LD    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/chrom2L_ld_window_avg.txt' delimiter = ',' MISSOVER
 DSD lrecl=32767 ;
       * informat size $4. ;
        informat chrom2l $6. ;
        informat bin best32. ;
        informat mess $6.   ;
        informat ld best32. ;
        *format size $4. ;
        format chrom2l $6. ;
        format bin best12. ;
        format mess $6.     ;
        format ld best12. ;
     input
                               chrom2l $
                 bin
                 mess $
                 ld
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


  data WORK.TD    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/all_chrom2L_10kb_window_tsd.txt' delimiter='09'x
 MISSOVER DSD lrecl=32767 firstobs=1 ;
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


data td1;
  set td;
  bin = bin_start/10000;
  run;

proc sort data=ld;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td;
  merge td1 (in=in1) ld (in=in2);
  by bin;
  run;

proc export data=ld_td
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chr2l_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;


/* Chromosome X */

 data WORK.TD_x    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/vcf_split_td_ld/all_chromX_tajimaD_10kb.Tajima.D' delimiter='09'x
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
     infile '/home/fnew/dsim/vcf_split_td_ld/chromX_ld_window_avg.txt' delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat dsim $4. ;
        informat window_size $4. ;
        informat chrom $6. ;
        informat bin best32. ;
        informat R2 best32. ;
        format dsim $4. ;
        format window_size $4. ;
        format chrom $6. ;
        format bin best12. ;
        format R2 best12. ;
     input
                 dsim $
                 window_size $
                 chrom $
                 bin
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
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chrX_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;




/* Chromosome 3L */

 data WORK.TD_3L    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/vcf_split_td_ld/all_chrom3L_tajimaD_10kb.Tajima.D' delimiter='09'x
  MISSOVER DSD lrecl=32767 ;
         informat chrom $2. ;
         informat bin_start best32. ;
         informat n_snps best32. ;
         informat TajimaD best32. ;
         format chrom $2. ;
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

 data WORK.LD_3L    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/chrom3L_ld_window_avg.txt' delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat window_size $4. ;
        informat chrom $7. ;
        informat bin best32. ;
        informat R2 best32. ;
        format window_size $4. ;
        format chrom $7. ;
        format bin best12. ;
        format R2 best12. ;
     input
                 window_size $
                 chrom $
                 bin
                 R2
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data td1;
  set td_3L;
  bin = bin_start/10000;
  run;

proc sort data=ld_3L;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td_3L;
  merge td1 (in=in1) ld_3L (in=in2);
  by bin;
  run;

proc export data=ld_td_3L
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chr3L_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;



/* Chromosome 3R */

 data WORK.TD_3R    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/vcf_split_td_ld/all_chrom3R_tajimaD_10kb.Tajima.D' delimiter='09'x
  MISSOVER DSD lrecl=32767 ;
         informat chrom $2. ;
         informat bin_start best32. ;
         informat n_snps best32. ;
         informat TajimaD best32. ;
         format chrom $2. ;
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

 data WORK.LD_3R    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/chrom3R_ld_window_avg.txt' delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat window_size $4. ;
        informat chrom $7. ;
        informat bin best32. ;
        informat R2 best32. ;
        format window_size $4. ;
        format chrom $7. ;
        format bin best12. ;
        format R2 best12. ;
     input
                 window_size $
                 chrom $
                 bin
                 R2
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data td1;
  set td_3R;
  bin = bin_start/10000;
  run;

proc sort data=ld_3R;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td_3R;
  merge td1 (in=in1) ld_3R (in=in2);
  by bin;
  run;

proc export data=ld_td_3R
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chr3R_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;


/* Chromosome 2R */

 data WORK.TD_2R    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/vcf_split_td_ld/all_chrom2R_tajimaD_10kb.Tajima.D' delimiter='09'x
  MISSOVER DSD lrecl=32767 ;
         informat chrom $2. ;
         informat bin_start best32. ;
         informat n_snps best32. ;
         informat TajimaD best32. ;
         format chrom $2. ;
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

 data WORK.LD_2R    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/vcf_split_td_ld/chrom2R_ld_window_avg.txt' delimiter = ','
 MISSOVER DSD lrecl=32767 ;
        informat window_size $4. ;
        informat chrom $7. ;
        informat bin best32. ;
        informat R2 best32. ;
        format window_size $4. ;
        format chrom $7. ;
        format bin best12. ;
        format R2 best12. ;
     input
                 window_size $
                 chrom $
                 bin
                 R2
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data td1;
  set td_2R;
  bin = bin_start/10000;
  run;

proc sort data=ld_2R;
  by bin;
  run;

proc sort data=td1;
  by bin;
  run;

data ld_td_2R;
  merge td1 (in=in1) ld_2R (in=in2);
  by bin;
  run;

proc export data=ld_td_2R
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/ld_and_D/chr2R_ld_and_td.txt"
    DBMS=TAB REPLACE;
    run;


