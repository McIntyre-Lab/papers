libname dsim "/home/fnew/dsim";

 data DSIM.allele_freqs    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/dmel/dmel_allele_freqs.frq' delimiter='09'x MISSOVER DSD lrecl=32767
 firstobs=2 ;
         informat CHROM $25. ;
         informat POS best32. ;
         informat N_ALLELES best32. ;
         informat N_CHR best32. ;
         informat allele_freq $15. ;
         informat allele_freq_2 $15. ;
         informat allele_freq_3 $15. ;
         informat allele_freq_4 $15.;
         format CHROM $25. ;
         format POS best12. ;
         format N_ALLELES best12. ;
         format N_CHR best12. ;
         format allele_freq $15. ;
         format allele_freq_2 $15.;
         format allele_freq_3 $15.;
         format allele_freq_4 $15.;
      input
                  CHROM $
                  POS
                  N_ALLELES
                  N_CHR
                  allele_freq $
                  allele_freq_2 $
                  allele_freq_3 $
                  allele_freq_4 $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


proc freq data=dsim.allele_freqs;
  tables n_alleles ; 
  run; *mono-, bi-, tri-, and tetra-allelic loci;
/* bi- 94.66%
   tri- 5.05%
   tet- 0.29%
   */


/* Import the frequencies by exon. */

   data DSIM.allele_freqs_exon    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/dsim_exons_allele_freq.frq' delimiter='09'x MISSOVER DSD lrecl=32767
 firstobs=2 ;
         informat CHROM $25. ;
         informat POS best32. ; 
         informat N_ALLELES best32. ;
         informat N_CHR best32. ;
         informat allele_freq $15. ;
         informat allele_freq_2 $15. ;
         informat allele_freq_3 $15. ;
         informat allele_freq_4 $15.;
         format CHROM $25. ;
         format POS best12. ;
         format N_ALLELES best12. ;
         format N_CHR best12. ;
         format allele_freq $15. ;
         format allele_freq_2 $15.;
         format allele_freq_3 $15.;
         format allele_freq_4 $15.;
      input
                  CHROM $
                  POS
                  N_ALLELES
                  N_CHR
                  allele_freq $
                  allele_freq_2 $
                  allele_freq_3 $
                  allele_freq_4 $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;
      
      
proc freq data=dsim.allele_freqs_exon;
  tables n_alleles ; 
  run;


  
/* Import the frequencies by intron. */

   data DSIM.allele_freqs_intron    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/dsim/dsim_introns_allele_freq.frq' delimiter='09'x MISSOVER DSD lrecl=32767
 firstobs=2 ;
         informat CHROM $25. ;
         informat POS best32. ; 
         informat N_ALLELES best32. ;
         informat N_CHR best32. ;
         informat allele_freq $15. ;
         informat allele_freq_2 $15. ;
         informat allele_freq_3 $15. ;
         informat allele_freq_4 $15.;
         format CHROM $25. ;
         format POS best12. ;
         format N_ALLELES best12. ;
         format N_CHR best12. ;
         format allele_freq $15. ;
         format allele_freq_2 $15.;
         format allele_freq_3 $15.;
         format allele_freq_4 $15.;
      input
                  CHROM $
                  POS
                  N_ALLELES
                  N_CHR
                  allele_freq $
                  allele_freq_2 $
                  allele_freq_3 $
                  allele_freq_4 $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;
      
      
proc freq data=dsim.allele_freqs_intron;
  tables n_alleles ; 
  run;
