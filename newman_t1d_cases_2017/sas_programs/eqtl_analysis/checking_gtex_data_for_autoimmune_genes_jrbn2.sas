ods listing; ods html close;
libname eqtl '!PATCON/eqtl_analysis/sas_data';
libname av '!PATCON/useful_human_data/aceview_hg19/sas_data';


/* 1. Get set of genes from the 9 GTEx pilot tissues that were tested in thier eqtl study
   2. flag gene*snp pairs with at least one significant eqtl
   3. Count ## shared between studies from all COMMON gene*snp pairs
	- this becomes the total %% of shared eqtl
   4. Convert the list of genes we testsed as eQTL from Aceview IDs to Ensembl IDs
   5. Subset GTEx results for genes we are interested in
   6. Calculate proportion of shared eqtl from this subset of genes
*/

/* Import multi-tissue eqtl results */
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

/*
    flag_multi_eqtl    Frequency     Percent     Frequency      Percent
    --------------------------------------------------------------------
                  0     4842543       45.92       4842543        45.92
                  1     5703041       54.08      10545584       100.00


                                               Cumulative    Cumulative
     flag_any_eqtl    Frequency     Percent     Frequency      Percent
     ------------------------------------------------------------------
                 1    10545584      100.00      10545584       100.00


                                               Cumulative    Cumulative
     flag_all_eqtl    Frequency     Percent     Frequency      Percent
     ------------------------------------------------------------------
                 0     2349106       22.28       2349106        22.28
                 1     8196478       77.72      10545584       100.00

*/


/* Import Ensembl2Entrez, then Ensembl2Aceview*/

proc import datafile="/home/jrbnewman/Downloads/gtex_data/ensembl2entrez_id.txt"
    out=ens2entrez dbms=tab replace; guessingrows=67043;
run;

data ens2entrez2;
  set ens2entrez;
  rename gene_stable_id=ens_id entrezgene_id=entrez_id;
run;

data av2entrez;
   set av.aceview2entrez;
   if entrez_id=. then delete;
run;

data eqtl_genes;
  set eqtl.snp2gene_index;
  keep gene_id;
run;

proc sort data=eqtl_genes nodup;
   by gene_id;
proc sort data=av2entrez nodup;
   by gene_id;
run;

data eqtl_av2entrez;
  merge eqtl_genes (in=in1) av2entrez (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=ens2entrez2 nodup;
   by entrez_id;
proc sort data=eqtl_av2entrez nodup;
  by entrez_id;
run;

data ens2av;
  merge ens2entrez2 (in=in1) eqtl_av2entrez (in=in2);
  by entrez_id;
  if in1 and in2;
run;

/* Subset eqtl */

data ens2av2;
  set ens2av;
  keep ens_id;
run;

data flag_multi2;
  length ens_id $18.;
  set flag_multi;
  ens_id=scan(gene,1,".");
run;

proc sort data=flag_multi2;
  by ens_id;
proc sort data=ens2av2 nodup;
  by ens_id;
run;

data flag_multi_immuno;
  merge flag_multi2 (in=in1) ens2av2 (in=in2);
  by ens_id;
  if in1 and in2;
run;



proc freq data=flag_multi_immuno;
  tables flag_multi_eqtl flag_any_eqtl flag_all_eqtl;
run;

/*
                                              Cumulative    Cumulative
  flag_multi_eqtl    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0      321384       54.86        321384        54.86
                1      264455       45.14        585839       100.00


                                             Cumulative    Cumulative
   flag_any_eqtl    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
               1      585839      100.00        585839       100.00


                                             Cumulative    Cumulative
   flag_all_eqtl    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
               0      126740       21.63        126740        21.63
               1      459099       78.37        585839       100.00


*/

data eqtl_genes2 eqtl_genes_nomult;
   length gene_id2 $36.;
   set eqtl.results_summary_table_w_means_v2;
   if flag_multigene=0 then do;
       gene_id2=gene_id;
       output eqtl_genes_nomult;
       end;
   else do;
      do i=1 by 1 while(scan(gene_id,i,"|") ^= "" );
          gene_id2=scan(gene_id,i,"|");
          output  eqtl_genes2;
          end;
     end;
   keep gene_id2;
  rename gene_id2=gene_id;
run;

proc sort data=eqtl_genes2 nodup;
  by gene_id;
proc sort data=eqtl_genes_nomult nodup;
  by gene_id;
run;


data eqtl_av2entrez;
  merge eqtl_genes_nomult (in=in1) av2entrez (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=ens2entrez2 nodup;
   by entrez_id;
proc sort data=eqtl_av2entrez nodup;
  by entrez_id;
run;

data ens2av;
  merge ens2entrez2 (in=in1) eqtl_av2entrez (in=in2);
  by entrez_id;
  if in1 and in2;
run;

/* Subset eqtl */

data ens2av2;
  set ens2av;
  keep ens_id;
run;

data flag_multi2;
  length ens_id $18.;
  set flag_multi;
  ens_id=scan(gene,1,".");
run;

proc sort data=flag_multi2;
  by ens_id;
proc sort data=ens2av2 nodup;
  by ens_id;
run;

data flag_multi_immuno;
  merge flag_multi2 (in=in1) ens2av2 (in=in2);
  by ens_id;
  if in1 and in2;
run;



proc freq data=flag_multi_immuno;
  tables flag_multi_eqtl flag_any_eqtl flag_all_eqtl;
run;

/*
SELECTED GENES

                                              Cumulative    Cumulative
  flag_multi_eqtl    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0      321384       54.86        321384        54.86
                1      264455       45.14        585839       100.00


                                             Cumulative    Cumulative
   flag_any_eqtl    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
               1      585839      100.00        585839       100.00


                                             Cumulative    Cumulative
   flag_all_eqtl    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
               0      126740       21.63        126740        21.63
               1      459099       78.37        585839       100.00


TESTED GENES


                                               Cumulative    Cumulative
   flag_multi_eqtl    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------------
                 0      202435       53.55        202435        53.55
                 1      175575       46.45        378010       100.00


                                              Cumulative    Cumulative
    flag_any_eqtl    Frequency     Percent     Frequency      Percent
    ------------------------------------------------------------------
                1      378010      100.00        378010       100.00


                                              Cumulative    Cumulative
    flag_all_eqtl    Frequency     Percent     Frequency      Percent
    ------------------------------------------------------------------
                0       84398       22.33         84398        22.33
                1      293612       77.67        378010       100.00


TESTED GENES (NO MULTI)


                                               Cumulative    Cumulative
   flag_multi_eqtl    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------------
                 0      318709       54.87        318709        54.87
                 1      262182       45.13        580891       100.00


                                              Cumulative    Cumulative
    flag_any_eqtl    Frequency     Percent     Frequency      Percent
    ------------------------------------------------------------------
                1      580891      100.00        580891       100.00


                                              Cumulative    Cumulative
    flag_all_eqtl    Frequency     Percent     Frequency      Percent
    ------------------------------------------------------------------
                0      125194       21.55        125194        21.55
                1      455697       78.45        580891       100.00



*/
