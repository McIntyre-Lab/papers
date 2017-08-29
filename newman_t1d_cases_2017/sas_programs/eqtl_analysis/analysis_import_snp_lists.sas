/* Import SNP lists */

/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Import genotyped SNP list */

    data WORK.GENOTYPED_SNPS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/home/jrbnewman/concannon/eqtl_analysis/original_data/eurfam_03202014.bim'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=1 ;
       informat chr $4. ;
       informat snp_id $15. ;
       informat centiMorgan best32. ;
       informat pos best32. ;
       informat minor_allele $1. ;
       informat major_allele $1. ;
       format chr $4. ;
       format snp_id $15.;
       format centiMorgan best12. ;
       format pos best12. ;
       format minor_allele $1. ;
       format major_allele $1. ;
    input
                chr
                snp_id $
                centiMorgan
                pos
                minor_allele $
                major_allele $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Import credible SNP table */

     data WORK.ONENGUT_SUPPTABLE1    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/jrbnewman/concannon/eqtl_analysis/snp_info/ic_eur_credible_snplist.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat chr $5. ;
        informat name $25. ;
        informat position best32. ;
        informat name_rs $11. ;
        informat snp_id $11. ;
        informat snp_type $8. ;
        informat credible_type $11. ;
        informat LZ_annotation $20. ;
        informat index_snp $16. ;
        informat index_snp_rs $11. ;
        informat index_snp_position best32. ;
        informat credible_snp $25. ;
        informat cred_snp_rs $11. ;
        informat credible_snp_pos best32. ;
        informat pp best32. ;
        informat r2 best32. ;
        informat MAF $11. ;
        informat func $10. ;
        informat genes $12. ;
        informat polyphen $5. ;
     informat sift $4. ;
     informat aa $9. ;
     informat enrich_enhancer_overlap best32. ;
     informat non_enrich_enhancer_overlap best32. ;
     format chr $5. ;
     format name $25. ;
     format position best12. ;
     format name_rs $11. ;
     format snp_id $11. ;
     format snp_type $8. ;
     format credible_type $11. ;
     format LZ_annotation $20. ;
     format index_snp $16. ;
     format index_snp_rs $11. ;
     format index_snp_position best12. ;
     format credible_snp $25. ;
     format cred_snp_rs $11. ;
     format credible_snp_pos best12. ;
     format pp best12. ;
     format r2 best12. ;
     format MAF $11. ;
     format func $10. ;
     format genes $12. ;
     format polyphen $5. ;
     format sift $4. ;
     format aa $9. ;
     format enrich_enhancer_overlap best12. ;
     format non_enrich_enhancer_overlap best12. ;
  input
               chr $
               name $
               position
               name_rs $
               snp_id $
               snp_type $
               credible_type $
               LZ_annotation $
               index_snp $
               index_snp_rs $
               index_snp_position
               credible_snp $
               cred_snp_rs $
               credible_snp_pos
               pp
               r2
               MAF $
               func $
               genes $
               polyphen $
               sift $
               aa $
               enrich_enhancer_overlap
               non_enrich_enhancer_overlap
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;


/* Make lists permenant for now */

data eqtl.genotyped_snps;
   set genotyped_snps;
run;

data eqtl.onengut_supptable1;
   set onengut_supptable1;
run;

/* Make working SNP lists for checking eQTLs */

data eqtl.supptable1_snp_list;
   set eqtl.onengut_supptable1;
   keep chr position snp_id index_snp_rs cred_snp_rs r2;
run;

