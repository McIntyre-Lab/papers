/* Check the overlap between T1D-association SNPs (and the statistically indistinguishable SNPs) and T1D splicing cis-EQTL */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/*
TO DO:
1. Import the "credible" SNP lists
2. Make a list of "tested SNPs"-to-"linked ImmunoChip SNPs"
3. Check overlap
*/


/* Import lists of credible SNPs from Onengut-Gumuscu et al, 2015 */
   These are excel formatted tables */

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


/* Make this table permenant for now */

data eqtl.onengut_supptable1;
   set onengut_supptable1;
run;

/* Make working SNP lists for checking eQTLs */


data eqtl.supptable1_snp_list;
   set eqtl.onengut_supptable1;
   keep chr position name_rs snp_id index_snp_rs cred_snp_rs r2;
run;



