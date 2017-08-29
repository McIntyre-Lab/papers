
libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Import merged VCF file */
/*
proc import datafile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/merged_vcf.tsv'
    out=snp_data_unformatted
    dbms=tab
    replace;
    guessingrows=12800;
    getnames=no;
run;

*/

/* Double check from the VCF (before extracting by gene) that these are the column headings as these will correspond to the variable names here: (first variable is GENE_ID!!)

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  M041_M041       M020_M020       M053_M053
       M007_M007       M027_M027       M018_M018       M005_M005       M006_M006       M016_M016       M070_M070
       M042_M042       M082_M082       M035_M035       M040_M040       M073_M073       M008_M008       M064_M064
       M033_M033       M052_M052       M066_M066       M036_M036       M065_M065       M063_M063       M062_M062
       M028_M028       M022_M022       M046_M046       M083_M083       M031_M031       M034_M034       M078_M078
       M039_M039       M074_M074       M058_M058       M001_M001       M080_M080       M072_M072       M055_M055
       M056_M056       M019_M019       M049_M049       M038_M038       M023_M023       M068_M068       M009_M009
       M002_M002       M050_M050       M021_M021       M012_M012       M077_M077       M048_M048       M026_M026
       M010_M010       M067_M067       M004_M004       M051_M051       M037_M037       M057_M057       M059_M059
       M081_M081       M017_M017       M014_M014       M043_M043       M025_M025       M060_M060       M029_M029
       M075_M075       M030_M030       M071_M071       M047_M047       M061_M061       M015_M015       M076_M076
       M069_M069       M003_M003       M054_M054       M013_M013       M011_M011       M044_M044       M024_M024
       M045_M045       M079_M079       M032_M032

*/

  /**********************************************************************
  *   PRODUCT:   SAS
  *   VERSION:   9.4
  *   CREATOR:   External File Interface
  *   DATE:      05NOV15
  *   DESC:      Generated SAS Datastep Code
  *   TEMPLATE SOURCE:  (None Specified.)
  ***********************************************************************/
     data WORK.SNP_DATA_UNFORMATTED    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/merged_vcf.tsv' delimiter='09'x MISSOVER DSD
 lrecl=32767 ;
 informat gene_id $33. ;   informat chr $4. ;      informat pos best32. ;    informat snp_id $15. ;     informat ref_allele $1. ;
 informat alt_allele $1. ;       informat qual best32. ; informat filter best32. ;   informat info best32. ;  informat format $2. ;
 informat M041 $3. ;      informat M020 $3. ;    informat M053 $3. ;      informat M007 $3. ;     informat M027 $3. ;
 informat M018 $3. ;      informat M005 $3. ;    informat M006 $3. ;      informat M016 $3. ;     informat M070 $3. ;
 informat M042 $3. ;      informat M082 $3. ;    informat M035 $3. ;      informat M040 $3. ;     informat M073 $3. ;
 informat M008 $3. ;      informat M064 $3. ;    informat M033 $3. ;      informat M052 $3. ;     informat M066 $3. ;
 informat M036 $3. ;      informat M065 $3. ;    informat M063 $3. ;      informat M062 $3. ;     informat M028 $3. ;
 informat M022 $3. ;      informat M046 $3. ;    informat M083 $3. ;      informat M031 $3. ;     informat M034 $3. ;
 informat M078 $3. ;      informat M039 $3. ;    informat M074 $3. ;      informat M058 $3. ;     informat M001 $3. ;
 informat M080 $3. ;      informat M072 $3. ;    informat M055 $3. ;      informat M056 $3. ;     informat M019 $3. ;
 informat M049 $3. ;      informat M038 $3. ;    informat M023 $3. ;      informat M068 $3. ;     informat M009 $3. ;
 informat M002 $3. ;      informat M050 $3. ;    informat M021 $3. ;      informat M012 $3. ;     informat M077 $3. ;
 informat M048 $3. ;      informat M026 $3. ;    informat M010 $3. ;      informat M067 $3. ;     informat M004 $3. ;
 informat M051 $3. ;      informat M037 $3. ;    informat M057 $3. ;      informat M059 $3. ;     informat M081 $3. ;
 informat M017 $3. ;      informat M014 $3. ;    informat M043 $3. ;      informat M025 $3. ;     informat M060 $3. ;
 informat M029 $3. ;      informat M075 $3. ;    informat M030 $3. ;      informat M071 $3. ;     informat M047 $3. ;
 informat M061 $3. ;      informat M015 $3. ;    informat M076 $3. ;      informat M069 $3. ;     informat M003 $3. ;
 informat M054 $3. ;      informat M013 $3. ;    informat M011 $3. ;      informat M044 $3. ;     informat M024 $3. ;
 informat M045 $3. ;      informat M079 $3. ;    informat M032 $3. ;
  format gene_id $29. ;   format chr $4. ;      format pos best12. ;    format snp_id $15. ;     format ref_allele $1. ;
 format alt_allele $1. ;       format qual best12. ; format filter best12. ;   format info best12. ;  format format $2. ;
 format M041 $3. ;      format M020 $3. ;    format M053 $3. ;      format M007 $3. ;     format M027 $3. ;
 format M018 $3. ;      format M005 $3. ;    format M006 $3. ;      format M016 $3. ;     format M070 $3. ;
 format M042 $3. ;      format M082 $3. ;    format M035 $3. ;      format M040 $3. ;     format M073 $3. ;
 format M008 $3. ;      format M064 $3. ;    format M033 $3. ;      format M052 $3. ;     format M066 $3. ;
 format M036 $3. ;      format M065 $3. ;    format M063 $3. ;      format M062 $3. ;     format M028 $3. ;
 format M022 $3. ;      format M046 $3. ;    format M083 $3. ;      format M031 $3. ;     format M034 $3. ;
 format M078 $3. ;      format M039 $3. ;    format M074 $3. ;      format M058 $3. ;     format M001 $3. ;
 format M080 $3. ;      format M072 $3. ;    format M055 $3. ;      format M056 $3. ;     format M019 $3. ;
 format M049 $3. ;      format M038 $3. ;    format M023 $3. ;      format M068 $3. ;     format M009 $3. ;
 format M002 $3. ;      format M050 $3. ;    format M021 $3. ;      format M012 $3. ;     format M077 $3. ;
 format M048 $3. ;      format M026 $3. ;    format M010 $3. ;      format M067 $3. ;     format M004 $3. ;
 format M051 $3. ;      format M037 $3. ;    format M057 $3. ;      format M059 $3. ;     format M081 $3. ;
 format M017 $3. ;      format M014 $3. ;    format M043 $3. ;      format M025 $3. ;     format M060 $3. ;
 format M029 $3. ;      format M075 $3. ;    format M030 $3. ;      format M071 $3. ;     format M047 $3. ;
 format M061 $3. ;      format M015 $3. ;    format M076 $3. ;      format M069 $3. ;     format M003 $3. ;
 format M054 $3. ;      format M013 $3. ;    format M011 $3. ;      format M044 $3. ;     format M024 $3. ;
 format M045 $3. ;      format M079 $3. ;    format M032 $3. ;
   input    gene_id $     chr $        pos       snp_id $       ref_allele $ 
  alt_allele $         qual    filter      info  format $ 
  M041 $        M020 $      M053 $        M007 $       M027 $ 
  M018 $        M005 $      M006 $        M016 $       M070 $ 
  M042 $        M082 $      M035 $        M040 $       M073 $ 
  M008 $        M064 $      M033 $        M052 $       M066 $ 
  M036 $        M065 $      M063 $        M062 $       M028 $ 
  M022 $        M046 $      M083 $        M031 $       M034 $ 
  M078 $        M039 $      M074 $        M058 $       M001 $ 
  M080 $        M072 $      M055 $        M056 $       M019 $ 
  M049 $        M038 $      M023 $        M068 $       M009 $ 
  M002 $        M050 $      M021 $        M012 $       M077 $ 
  M048 $        M026 $      M010 $        M067 $       M004 $ 
  M051 $        M037 $      M057 $        M059 $       M081 $ 
  M017 $        M014 $      M043 $        M025 $       M060 $ 
  M029 $        M075 $      M030 $        M071 $       M047 $ 
  M061 $        M015 $      M076 $        M069 $       M003 $ 
  M054 $        M013 $      M011 $        M044 $       M024 $ 
  M045 $        M079 $      M032 $ 
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;


/* untranspose data */
    
data snp_data_unformatted2;
   set snp_data_unformatted;
run;

proc sort data=snp_data_unformatted2;
   by gene_id snp_id;
run;

proc transpose
    Name=subject_id
    data=snp_data_unformatted2
    out=snp_data_tposed(rename=(col1=genotype));
    by gene_id snp_id;
    var M041 M020 M053 M007  M027 M018 M005 M006 M016 M070 M042 M082 M035 M040 M073 M008 M064 M033 M052 M066 M036 
        M065 M063 M062 M028 M022 M046 M083 M031 M034 M078 M039 M074 M058 M001 M080 M072 M055 M056 M019 M049 M038 
        M023 M068 M009 M002 M050 M021 M012 M077 M048 M026 M010 M067 M004 M051 M037 M057 M059 M081 M017 M014 M043 
       M025 M060 M029 M075 M030 M071 M047 M061 M015 M076 M069 M003 M054 M013 M011 M044 M024 M045 M079 M032;
run;


/* Change genotypes to 0,1,2 */

data snp_data_tposed2;
set snp_data_tposed;
   if genotype='./.' then genotype012=.;
   else if genotype='0/0' then genotype012=0;
   else if genotype='0/1' then genotype012=1;
   else if genotype='1/1' then genotype012=2;
   else genotype012=.;
   drop genotype;
   rename genotype012=genotype;
run;

/* Get SNP info */

data snp_info;
   set snp_data_unformatted2;
   keep gene_id chr pos	snp_id ref_allele alt_allele;
run;

proc sort data=snp_info nodup;
   by gene_id snp_id;
run;

proc sort data=snp_data_tposed2;
   by gene_id snp_id;
run;

data snp_data_w_info oops1 oops2;
   merge snp_info (in=in1) snp_data_tposed2 (in=in2);
   by gene_id snp_id;
   if in1 and in2 then output snp_data_w_info;
   else if in1 then output oops1;
   else output oops2;
run;

/* Make permenant */

data eqtl.snp_data_w_info;
   set snp_data_w_info;
run;

