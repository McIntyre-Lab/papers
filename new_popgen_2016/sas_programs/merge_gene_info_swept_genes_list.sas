/* I have created a list of genes that overlap the top50 H12 peaks from each
 * chromosome. I want to now merge in gene info and GO terms */

libname sim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";


  data WORK.Swept_genes    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/top_peaks_swept_genes_allchrom.txt' delimiter='09'x MISSOVER DSD lrecl=32767 ;
         informat gene_chrom $2. ;
         informat gene_start best32. ;
         informat gene_end best32. ;
         informat fbgn $11. ;
         informat VAR5 best32. ;
         informat strand $1. ;
         informat chrom $2. ;
         informat peak_start best32. ;
         informat peak_end best32. ;
         informat H12 best32. ;
         format gene_chrom $2. ;
         format gene_start best12. ;
         format gene_end best12. ;
         format fbgn $11. ;
         format VAR5 best12. ;
         format strand $1. ;
         format chrom $2. ;
         format peak_start best12. ;
         format peak_end best12. ;
         format H12 best12. ;
      input
                  gene_chrom $
                  gene_start
                  gene_end
                  fbgn $
                  VAR5
                  strand $
                  chrom $
                  peak_start
                  peak_end
                  H12
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

*Import the output from flymine;
  data WORK.Sim2mel    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/Downloads/swept_genes_dsim.tsv' delimiter='09'x MISSOVER DSD
 lrecl=32767 firstobs=1 ;
        informat gene $7. ;
        informat sim_symbol $19. ;
        informat FBgn $11. ;
        informat species $19. ;
        format gene $7. ;
        format sim_symbol $19. ;
        format FBgn $11. ;
        format species $19. ;
     input
                 gene $
                 sim_symbol $
                 FBgn $
                 species $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;




*merge them together;
proc sort data=swept_genes; by fbgn; run;
proc sort data=sim2mel; by fbgn; run;

data sim_swept_genes;
  merge swept_genes (in=in1) sim2mel (in=in2);
  by fbgn;
  if in1 and in2;
  run;

proc export data=sim_swept_genes
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/H12/swept_genes_top_peaks.csv"
    DBMS=CSV REPLACE;
    run;
