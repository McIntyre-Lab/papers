/* Check the overlap between T1D-association SNPs (and the statistically indistinguishable SNPs) and T1D splicing cis-EQTL */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';


/* Import tested SNP lists and SNPs tagged by them
   Will use this to convert credible SNPs to tested SNPs */

/* PLINK output is fixed with outputs so I will need to import them differently */

* Set filenames;
filename pruned '/home/jrbnewman/concannon/eqtl_analysis/original_data/pruned_list_ld_all.ld';
filename tagged '/home/jrbnewman/concannon/eqtl_analysis/original_data/tagging_snps.tags.list';
filename kept '/home/jrbnewman/concannon/eqtl_analysis/original_data/pruned_snps.prune.in';

/* Import LD data */
data snp_ld_data;
infile pruned truncover;
input snp_id_1 $ 21-35 snp_id_2 $ 57-71 r2 $ 73-85;
run;

/* Import tagging SNPs */
data tagging_snps;
infile tagged truncover;
input snp_id $ 1-15 tagged_snps $ 69-3000;
run;

/* Import tagging SNPs */
    data WORK.TAGGING_SNPS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/home/jrbnewman/concannon/eqtl_analysis/original_data/tagging_snp_list.csv'
delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat snp_id $15. ;
       informat tagged_snps $1618. ;
       format snp_id $15. ;
       format tagged_snps $1618. ;
    input
                snp_id $
                tagged_snps $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Import list of kept SNPs */
data kept_snps;
infile kept truncover;
input snp_id $ 1-15;
run;


/* Drop first row: header row from PLINK output
   Strip leading and trailing spaces */



data snp_ld_data_2;
  set snp_ld_data;
  if _n_ = 1 then delete;
  snp_a_str=strip(put(snp_id_1, 15.));
  snp_b_str=strip(put(snp_id_2, 15.));
  r2_num=r2*1;
  keep snp_a_str snp_b_str r2_num;
  rename snp_a_str=snp_id_1;
  rename snp_b_str=snp_id_2;
  rename r2_num=r2;
run;

data kept_snps_2;
   set kept_snps;
   snp_id_str=strip(put(snp_id, 15.));
   keep snp_id_str;qq
   rename snp_id_str=snp_id;
run;

/* Make imported data permenant */

data eqtl.tagged_snp_list;
   set tagging_snps;
run;

data eqtl.filtered_snp_list;
   set kept_snps_2;
run;

data eqtl.snp_ld_data;
   set snp_ld_data_2;
run;
