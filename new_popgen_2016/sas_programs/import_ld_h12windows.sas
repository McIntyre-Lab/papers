/* Import the LD by h12 windows data for each chromosome */

 libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

  data WORK.CHR2L    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chr2L_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


 data WORK.CHR2R    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chr2R_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;     

 data WORK.CHR3R    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chr3R_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;    

data WORK.CHR3L    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chr3L_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;     


data WORK.CHR4    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chr4_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;     


data WORK.CHRX    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/LD/LD_in_h12_windows/chrX_ld_h12window_avg.csv' delimiter = ',' MISSOVER DSD lrecl=32767 ;
        informat Average $7. ;
        informat peak $4. ;
        informat center best32. ;
        informat LD best32. ;
        format Average $7. ;
        format peak $4. ;
        format center best12. ;
        format LD best12. ;
     input
                 Average $
                 peak $
                 center
                 LD
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;     

data dsim.chr2L_LD_h12windows;
  set chr2l;
  run;

data dsim.chr2R_LD_h12windows;
  set chr2r;
  run;

data dsim.chr3L_LD_h12windows;
  set chr3l;
  run;

data dsim.chr3R_LD_h12windows;
  set chr3R;
  run;

data dsim.chr4_LD_h12windows;
  set chr4;
  run;

data dsim.chrX_LD_h12windows;
  set chrx;
  run;


