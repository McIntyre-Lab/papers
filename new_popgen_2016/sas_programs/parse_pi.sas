/* Parse pi output for plotting */

libname dsim "/home/fnew/mclab/ethanol/Sim_Pop_Gen/sas_data";
libname fnew "/home/fnew/dsim";


     data WORK.pi    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/dsim_pi_1000bpwindow.windowed.pi' delimiter='09'x MISSOVER DSD
 lrecl=32767 firstobs=2 ;
         informat CHROM $20. ;
         informat BIN_START best32. ;
         informat BIN_END best32. ;
         informat N_VARIANTS best32. ;
         informat PI best32. ;
         format CHROM $20. ;
         format BIN_START best12. ;
         format BIN_END best12. ;
         format N_VARIANTS best12. ;
         format PI best12. ;
      input
                  CHROM $
                  BIN_START
                  BIN_END
                  N_VARIANTS
                  PI
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


*Separate by chromosomes and do X/A;

data chrx;
    set pi;
    if chrom="X";
    run;

data autosomes;
    set pi;
    if chrom ne "X";
    run;

data chr2r;
    set pi;
    if chrom="2R";
    run;
data chr2L;
    set pi;
    if chrom="2L";
    run;
data chr3R;
    set pi;
    if chrom="3R";
    run;
data chr3L;
    set pi;
    if chrom="3L";
    run;

data contigs;
    set pi;
    if chrom like "NODE";
    run;

