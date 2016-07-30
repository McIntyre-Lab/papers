libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

data dsim.Tsd_10kb    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile ‘!MCLAB/ethanol/Sim_Pop_Gen/output/tajima_d/filter_nolab_tajd_10kb.Tajima.D' delimiter='09'x
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


*Need only the significant windows, add a window end;
data signif;
  retain chrom BIN_START BIN_END ;
  set dsim.tsd_10kb;
  if TajimaD > 2.095 or TajimaD < -1.765;
  BIN_END=BIN_START+9999;
  run; *1080 obs;



proc export data=signif
    outfile=“!MCLAB/ethanol/Sim_Pop_Gen/output/tajima_d/significant_windows_10kb_tsd.txt"
    DBMS=TAB REPLACE;
    run;
