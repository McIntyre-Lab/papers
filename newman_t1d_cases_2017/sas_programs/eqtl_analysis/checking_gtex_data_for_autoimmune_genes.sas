ods listing; ods html close;
libname eqtl '!PATCON/eqtl_analysis/sas_data';


/* 1. Get set of genes from the 9 GTEx pilot tissues that were tested in thier eqtl study
   2. flag gene*snp pairs with at least one significant eqtl
   3. Count ## shared between studies from all COMMON gene*snp pairs
	- this becomes the total %% of shared eqtl
   4. Convert the list of genes we testsed as eQTL from Aceview IDs to Ensembl IDs
   5. Subset GTEx results for genes we are interested in
   6. Calculate proportion of shared eqtl from this subset of genes
*/

/* 1. Get set of genes from the 9 GTEx pilot tissues that were tested in thier eqtl study */

%macro importGenes(datain,tissueID);

proc import datafile="/home/jrbnewman/Downloads/gtex_data/&datain." out=&tissueID._expr
     dbms=tab replace; guessingrows=28344;
run;

data &tissueID._genes;
   set &tissueID._expr;
   keep Id;
run;

proc sort data=&tissueID._genes nodup;
   by Id;
run;

%mend;

%importGenes(Adipose_Subcutaneous.expr.txt,adipose);
%importGenes(Artery_Tibial.expr.txt,artery);
%importGenes(Heart_Left_Ventricle.expr.txt,heart);
%importGenes(Lung.expr.txt,lung);
%importGenes(Muscle_Skeletal.expr.txt,muscle);
%importGenes(Nerve_Tibial.expr.txt,nerve);
%importGenes(Skin_Sun_Exposed_Lower_leg.expr.txt,skin);
%importGenes(Thyroid.expr.txt,thyroid);
%importGenes(Whole_Blood.expr.txt,blood);

*merge all together

data gtex_genes;
  merge adipose_genes (in=in1) artery_genes(in=in2) heart_genes(in=in3) lung_genes(in=in4)
        muscle_genes(in=in5) nerve_genes(in=in6) skin_genes(in=in7) thyroid_genes(in=in8) blood_genes(in=in9);
  by id;
  if in1 then flag_adipose=1; else flag_adipose=0;
  if in2 then flag_artery=1; else flag_artery=0;
  if in3 then flag_heart=1; else flag_heart=0;
  if in4 then flag_lung=1; else flag_lung=0;
  if in5 then flag_muscle=1; else flag_muscle=0;
  if in6 then flag_nerve=1; else flag_nerve=0;
  if in7 then flag_skin=1; else flag_skin=0;
  if in8 then flag_thyroid=1; else flag_thyroid=0;
  if in9 then flag_blood=1; else flag_blood=0;
run;
proc freq data=gtex_genes noprint;
   tables flag_adipose*flag_artery*flag_heart*flag_lung*flag_muscle*flag_nerve*flag_skin*flag_thyroid*flag_blood / out=gene_count;
proc print data=gene_count;
run; *19039 genes in all 9 tissues;

/* 2. Flag any gene-snp pair that is signficant */


%macro importResults(datain,tissueID);

    data WORK.&tissueID._eqtl    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/home/jrbnewman/Downloads/gtex_data/&datain."
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat SNP $21. ;
       informat gene $18. ;
       informat t_stat best32. ;
       informat p_value best32. ;
       informat FDR best32. ;
       informat nom_thresh best32. ;
       informat min_p_ best32. ;
       informat EmpP best32. ;
       informat Ks best32. ;
       informat n best32. ;
       format SNP $21. ;
       format gene $18. ;
       format t_stat best12. ;
       format p_value best12. ;
       format FDR best12. ;
       format nom_thresh best12. ;
       format min_p_ best12. ;
       format EmpP best12. ;
       format Ks best12. ;
       format n best12. ;
    input
                SNP $
                gene $
                t_stat
                p_value
                FDR
                nom_thresh
                min_p_
                EmpP
                Ks
                n
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

data &tissueID._pairs;
   set &tissueID._eqtl;
   keep gene snp;
   rename gene=id;
run;

proc sort data=&tissueID._pairs nodup;
   by Id snp;
run;

%mend;

%importResults(Adipose_Subcutaneous.snp-gene-eQTLs.txt,adipose);
%importResults(Artery_Tibial.snp-gene-eQTLs.txt,artery);
%importResults(Heart_Left_Ventricle.snp-gene-eQTLs.txt,heart);
%importResults(Lung.snp-gene-eQTLs.txt,lung);
%importResults(Muscle_Skeletal.snp-gene-eQTLs.txt,muscle);
%importResults(Nerve_Tibial.snp-gene-eQTLs.txt,nerve);
%importResults(Skin_Sun_Exposed_Lower_leg.snp-gene-eQTLs.txt,skin);
%importResults(Thyroid.snp-gene-eQTLs.txt,thyroid);
%importResults(Whole_Blood.snp-gene-eQTLs.txt,blood);



data gtex_eqtl;
  merge adipose_pairs (in=in1) artery_pairs(in=in2) heart_pairs(in=in3) lung_pairs(in=in4)
        muscle_pairs(in=in5) nerve_pairs(in=in6) skin_pairs(in=in7) thyroid_pairs(in=in8) blood_pairs(in=in9);
  by id snp;
  if in1 then flag_adipose_eqtl=1; else flag_adipose_eqtl=0;
  if in2 then flag_artery_eqtl=1; else flag_artery_eqtl=0;
  if in3 then flag_heart_eqtl=1; else flag_heart_eqtl=0;
  if in4 then flag_lung_eqtl=1; else flag_lung_eqtl=0;
  if in5 then flag_muscle_eqtl=1; else flag_muscle_eqtl=0;
  if in6 then flag_nerve_eqtl=1; else flag_nerve_eqtl=0;
  if in7 then flag_skin_eqtl=1; else flag_skin_eqtl=0;
  if in8 then flag_thyroid_eqtl=1; else flag_thyroid_eqtl=0;
  if in9 then flag_blood_eqtl=1; else flag_blood_eqtl=0;
  sum_flags=flag_adipose_eqtl+flag_artery_eqtl+flag_heart_eqtl+flag_lung_eqtl+
          flag_muscle_eqtl+flag_nerve_eqtl+flag_skin_eqtl+flag_thyroid_eqtl+flag_blood_eqtl;
run;

proc freq data=gtex_eqtl noprint;
   tables flag_adipose_eqtl*flag_artery_eqtl*flag_heart_eqtl*flag_lung_eqtl*
          flag_muscle_eqtl*flag_nerve_eqtl*flag_skin_eqtl*flag_thyroid_eqtl*flag_blood_eqtl / out=eqtl_count;
proc print data=eqtl_count;
run;

proc freq data=gtex_eqtl;
   tables sum_flags;
run;

/* Merge eqtl with gene list */

proc sort data=gtex_eqtl;
  by id;
proc sort data=gtex_Genes;
  by id;
run;

data gtex_gene_w_eqtl;
  merge gtex_genes (in=in1) gtex_eqtl (in=in2);
  by id;
  if in1 and in2;
run;

/* Look at common genes only */

data gtex_gene_W_eqtl_common;
   set gtex_gene_W_eqtl;
   where flag_adipose=1 and 
         flag_artery=1 and 
         flag_heart=1 and 
         flag_lung=1 and 
         flag_muscle=1 and 
         flag_nerve=1 and 
         flag_skin=1 and 
         flag_thyroid=1 and 
         flag_blood=1 ;
run;

proc freq data=gtex_gene_W_eqtl_common;
   tables sum_flags;
run;

/* Try importing the multigene data instead */

    data WORK.ALL_SNPS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/Downloads/gtex_data/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_amea
n_com_genes_com_snps_all.txt' delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat gene $18. ;
       informat snp $21. ;
       informat Adipose best32. ;
       informat Artery best32. ;
       informat Blood best32. ;
       informat Heart best32. ;
       informat Lung best32. ;
       informat Muscle best32. ;
       informat Nerve best32. ;
       informat Skin best32. ;
       informat Thyroid best32. ;
       format gene $18. ;
       format snp $21. ;
       format Adipose best12. ;
       format Artery best12. ;
       format Blood best12. ;
       format Heart best12. ;
       format Lung best12. ;
       format Muscle best12. ;
       format Nerve best12. ;
       format Skin best12. ;
       format Thyroid best12. ;
    input
                gene $
                snp $
                Adipose
                Artery
                Blood
                Heart
                Lung
                Muscle
                Nerve
                Skin
                Thyroid
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


proc means data=all_snps noprint;
   var adipose artery blood heart lung muscle nerve skin thyroid;
   output out=check max=;
run;

data flag_check;
   set all_snps;
   if adipose=1 and artery=1 and blood=1 and heart=1 and lung=1 and muscle=1 and nerve=1 and skin=1 and thyroid=1
     then flag_all_one=1; else flag_all_one=0;
   if adipose=1 or artery=1 or blood=1 or heart=1 or lung=1 or muscle=1 or nerve=1 or skin=1 or thyroid=1
     then flag_any_one=1; else flag_any_one=0;
run;

proc freq data=flag_check;
  tables flag_all_one*flag_any_one flag_all_one flag_any_one;
run;


data flag_check;
   set all_snps;
   if adipose>0.8 and artery>0.8 and blood>0.8 and heart>0.8 and lung>0.8 and muscle>0.8 and nerve>0.8 and skin>0.8 and thyroid>0.8
     then flag_all_one=1; else flag_all_one=0;
   if adipose>0.8 or artery>0.8 or blood>0.8 or heart>0.8 or lung>0.8 or muscle>0.8 or nerve>0.8 or skin>0.8 or thyroid>0.8
     then flag_any_one=1; else flag_any_one=0;
run;

proc freq data=flag_check;
  tables flag_all_one*flag_any_one flag_all_one flag_any_one;
run;


proc import datafile='/home/jrbnewman/Downloads/gtex_data/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_amean_com_genes_com_snps_configs_all.txt'
  out=multi_eqtl dbms=tab replace;
run;

     data WORK.MULTI_EQTL    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/jrbnewman/Downloads/gtex_data/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_amean_com_genes_com_snps_configs_all.txt' delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat gene $18. ;
       informat snp $21. ;
        informat _000000001 best32. ;
        informat _000000010 best32. ;
        informat _000000100 best32. ;
        informat _000001000 best32. ;
        informat _000010000 best32. ;
        informat _000100000 best32. ;
        informat _001000000 best32. ;
        informat _010000000 best32. ;
        informat _100000000 best32. ;
        informat _111111111 best32. ;
        format gene $18. ;
        format snp $21. ;
        format _000000001 best12. ;
        format _000000010 best12. ;
        format _000000100 best12. ;
        format _000001000 best12. ;
        format _000010000 best12. ;
        format _000100000 best12. ;
        format _001000000 best12. ;
        format _010000000 best12. ;
        format _100000000 best12. ;
        format _111111111 best12. ;
   input
               gene $
               snp $
               _000000001
               _000000010
               _000000100
               _000001000
               _000010000
               _000100000
               _001000000
               _010000000
               _100000000
               _111111111
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

/* Flag multi-tissue eqtl */

data flag_multi;
  set multi_eqtl;
  if _111111111 < 0.05 then flag_multi_eqtl=1;
  else flag_multi_eqtl=0;
run;

data genes_on_all;
   set gtex_Genes;
   where flag_adipose=1 and 
         flag_artery=1 and 
         flag_heart=1 and 
         flag_lung=1 and 
         flag_muscle=1 and 
         flag_nerve=1 and 
         flag_skin=1 and 
         flag_thyroid=1 and 
         flag_blood=1 ;
   keep id;
   rename id=gene;
run;

proc sort data=flag_multi;
  by gene;
proc sort data=genes_on_all; 
  by gene;
run;

data genes_on_all_flag_multi;
  merge flag_multi (in=in1) genes_on_all (in=in2);
  by gene;
  if in1 and in2;
run;

proc freq data=genes_on_all_flag_multi;
  tables flag_multi_eqtl;
run;

proc means data=genes_on_all_flag_multi noprint;
   by gene;
   var flag_multi_eqtl;
   output out=flag_gene_multi_eqtl max=;
run;


proc freq data=flag_gene_multi_eqtl;
  tables flag_multi_eqtl;
run;

data flag_multi2;
  set flag_multi;
  if _000000001 < 0.05 
  or _000000010 < 0.05 
  or _000000100 < 0.05 
  or _000001000 < 0.05 
  or _000010000 < 0.05 
        or          _000100000 < 0.05 
        or          _001000000 < 0.05 
        or          _010000000 < 0.05 
        or          _100000000 < 0.05 
then flag_any_eqtl=1; else flag_any_eqtl=0;
run;


proc sort data=flag_multi2;
  by gene;
proc means data=flag_multi2 noprint;
   by gene;
   var flag_multi_eqtl flag_any_eqtl;
   output out=flag_gene_multi_eqtl3 max=;
run;


proc freq data=flag_gene_multi_eqtl3;
  tables flag_multi_eqtl*flag_any_eqtl flag_multi_eqtl flag_any_eqtl;
run;





     data WORK.MULTI_EQTL    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/jrbnewman/Downloads/gtex_data/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_unc_com_genes_com_snps_configs_all_sorted.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat gene $18. ;
       informat snp $21. ;
        informat _000000001 best32. ;
        informat _000000010 best32. ;
        informat _000000100 best32. ;
        informat _000001000 best32. ;
        informat _000010000 best32. ;
        informat _000100000 best32. ;
        informat _001000000 best32. ;
        informat _010000000 best32. ;
        informat _100000000 best32. ;
        informat _111111111 best32. ;
        format gene $18. ;
        format snp $21. ;
        format _000000001 best12. ;
        format _000000010 best12. ;
        format _000000100 best12. ;
        format _000001000 best12. ;
        format _000010000 best12. ;
        format _000100000 best12. ;
        format _001000000 best12. ;
        format _010000000 best12. ;
        format _100000000 best12. ;
        format _111111111 best12. ;
   input
               gene $
               snp $
               _000000001
               _000000010
               _000000100
               _000001000
               _000010000
               _000100000
               _001000000
               _010000000
               _100000000
               _111111111
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

/* Flag multi-tissue eqtl */

data flag_multi;
  set multi_eqtl;
  if _111111111 < 0.05 then flag_multi_eqtl=1;
  else flag_multi_eqtl=0;

  if _000000001 < 0.05 
  or _000000010 < 0.05 
  or _000000100 < 0.05 
  or _000001000 < 0.05 
  or _000010000 < 0.05 
        or          _000100000 < 0.05 
        or          _001000000 < 0.05 
        or          _010000000 < 0.05 
        or          _100000000 < 0.05 
then flag_any_eqtl=1; else flag_any_eqtl=0;

  if _000000001 < 0.05 
  and _000000010 < 0.05 
  and _000000100 < 0.05 
  and _000001000 < 0.05 
  and _000010000 < 0.05 
        and          _000100000 < 0.05 
        and          _001000000 < 0.05 
        and          _010000000 < 0.05 
        and          _100000000 < 0.05 
then flag_all_eqtl=1; else flag_all_eqtl=0;

run;


proc freq data=flag_multi;
  tables flag_multi_eqtl flag_any_eqtl flag_all_eqtl;
run;


proc sort data=flag_multi;
  by gene;
proc means data=flag_multi noprint;
   by gene;
   var flag_multi_eqtl flag_any_eqtl flag_all_eqtl;
   output out=flag_gene_multi_eqtl max=;
run;


proc freq data=flag_gene_multi_eqtl;
  tables flag_multi_eqtl flag_any_eqtl flag_all_eqtl;
run;
