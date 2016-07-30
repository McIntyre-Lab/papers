/* Not able to get theta from D and pi. I will try using SNP density. 

Miguel wrote a python script to get a1, so now i just need to divide s/a1

*/

libname dsim "/home/fnew/dsim/sas_data";


     data WORK.SNPDEN    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/snp_den/snpden_10kb_genome_filter_nolab_a1.tsv' delimiter='09'x
 MISSOVER DSD lrecl=32767 ;
        informat chrom $2. ;
        informat bin_start best32. ;
        informat n_snps best32. ;
        informat var_per_kb best32. ;
        informat a1 best32. ;
        format chrom $2. ;
        format bin_start best12. ;
        format n_snps best12. ;
        format var_per_kb best12. ;
        format a1 best12. ;
     input
                 chrom $
                 bin_start
                 n_snps
                 var_per_kb
                 a1
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data dsim.snpden;
  set snpden;
  run;

I:
data snpden2;
  set snpden;
  theta=n_snps/a1;
  run;

proc export data=snpden2
  outfile="/home/fnew/dsim/theta_from_snpden.txt"
  dbms=TAB REPLACE;
  run;

/* File is located: MCLAB/ethanol/Sim_Pop_Gen/output/snp_density/snpden_10kb_genome_filter_nolab_theta_4-20.tsv */