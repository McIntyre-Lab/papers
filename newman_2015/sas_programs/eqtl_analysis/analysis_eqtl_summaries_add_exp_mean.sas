/* Add expression means and ratios to eQTL results */

libname eqtl '/mnt/data/eqtls/sas_data2';

/* We want to annotate:
eqtl.eqtl_results_w_onengut_credible
eqtl.eqtl_results_summary_table */



/* Need to redo the expression means strings, as these seem to have not formatted completely correctly */


data genotype_means;
   set eqtl.mean_Coverage_per_genotype;
   format hmz1_mean_str $10.;
   format hmz1_sd_str $10.;
   format htz_mean_str $10.;
   format htz_sd_str $10.;
   format hmz2_mean_str $10.;
   format hmz2_sd_str $10.;
   hmz1_mean_str=strip(put(hmz1_mean, 10.3));
   hmz1_sd_str=strip(put(hmz1_sd, 10.3));
   htz_mean_str=strip(put(htz_mean, 10.3));
   htz_sd_str=strip(put(htz_sd, 10.3));
   hmz2_mean_str=strip(put(hmz2_mean, 10.3));
   hmz2_sd_str=strip(put(hmz2_sd, 10.3));
   keep feature_id snp_id cell_type hmz1_mean hmz1_sd htz_mean htz_sd hmz2_mean hmz2_sd
   hmz1_mean_str hmz1_sd_str htz_mean_str htz_sd_str hmz2_mean_str hmz2_sd_str
   htz_hmz1_ratio hmz2_hmz1_ratio ;
run;

/* Drop for duplicate observations */

proc sort data=genotype_means nodup;
   by feature_id snp_id cell_type;
run;


data genotype_means2;
   set genotype_means;
   length hmz1_expression $30.;
   length htz_expression $30.;
   length hmz2_expression $30.;
   if hmz1_mean=. then hmz1_expression='n.d.';
   else do;
       if hmz1_sd=. then hmz1_expression=strip(hmz1_mean_str);
       else hmz1_expression = cat(strip(hmz1_mean_str), ' ± ', strip(hmz1_sd_str));
       end;

   if htz_mean=. then htz_expression='n.d.';
   else do;
       if htz_sd=. then htz_expression=strip(htz_mean_str);
       else htz_expression = cat(strip(htz_mean_str), ' ± ', strip(htz_sd_str));
       end;

   if hmz2_mean=. then hmz2_expression='n.d.';
   else do;
       if hmz2_sd=. then hmz2_expression=strip(hmz2_mean_str);
       else hmz2_expression = cat(strip(hmz2_mean_str), ' ± ', strip(hmz2_sd_str));
       end;
   drop hmz1_mean hmz1_sd htz_mean htz_sd hmz2_mean hmz2_sd 
        hmz1_mean_str hmz1_sd_str htz_mean_str htz_sd_str hmz2_mean_str hmz2_sd_str; 
run;


/* Now I need to move each cell type to a side-by-side because of the way the eQTL summaries are formatted */

data cd4_means cd8_means cd19_means;
   set genotype_means2;
   if cell_type='CD4' then output cd4_means;
   if cell_type='CD8' then output cd8_means;
   if cell_type='CD19' then output cd19_means;
   drop cell_type;
run;

data cd4_means_2;
   set cd4_means;
   rename hmz1_expression=cd4_hmz1_expression
          htz_expression=cd4_htz_expression
          hmz2_expression=cd4_hmz2_expression
          htz_hmz1_ratio=cd4_hmz1_htz_hmz1_ratio
          hmz2_hmz1_ratio=cd4_hmz2_hmz1_ratio;
run;

data cd8_means_2;
   set cd8_means;
   rename hmz1_expression=cd8_hmz1_expression
          htz_expression=cd8_htz_expression
          hmz2_expression=cd8_hmz2_expression
          htz_hmz1_ratio=cd8_hmz1_htz_hmz1_ratio
          hmz2_hmz1_ratio=cd8_hmz2_hmz1_ratio;
run;

data cd19_means_2;
   set cd19_means;
   rename hmz1_expression=cd19_hmz1_expression
          htz_expression=cd19_htz_expression
          hmz2_expression=cd19_hmz2_expression
          htz_hmz1_ratio=cd19_hmz1_htz_hmz1_ratio
          hmz2_hmz1_ratio=cd19_hmz2_hmz1_ratio;
run;

proc sort data=cd4_means_2;
   by feature_id snp_id;
proc sort data=cd8_means_2;
   by feature_id snp_id;
proc sort data=cd19_means_2;
   by feature_id snp_id;
run;

data genotype_means_all;
   merge cd4_means_2 (in=in1) cd8_means_2 (in=in2) cd19_means_2 (in=in3);
   by feature_id snp_id;
   if not in1 then do;
          cd4_hmz1_expression='n.d.';
          cd4_htz_expression='n.d.';
          cd4_hmz2_expression='n.d.'; end;
   if not in2 then do;
          cd8_hmz1_expression='n.d.';
          cd8_htz_expression='n.d.';
          cd8_hmz2_expression='n.d.'; end;
   if not in3 then do;
          cd19_hmz1_expression='n.d.';
          cd19_htz_expression='n.d.';
          cd19_hmz2_expression='n.d.'; end;
run;


/* Merge genotype means into summary tables */

proc sort data=genotype_means_all;
   by feature_id snp_id;
proc sort data=eqtl.eqtl_results_summary_table;
   by feature_id snp_id;
proc sort data=eqtl.eqtl_results_w_onengut_credible;
   by feature_id snp_id;
run;


data eqtl_results_w_geno_means;
   merge eqtl.eqtl_results_summary_table (in=in1) genotype_means_all (in=in2);
   by feature_id snp_id;
   if in1 and in2;
run;

data eqtl_og_results_w_geno_means;
   merge eqtl.eqtl_results_w_onengut_credible (in=in1) genotype_means_all (in=in2);
   by feature_id snp_id;
   if in1 and in2;
run;



/* Make permenant */

data eqtl.results_summary_table_w_means;
   set eqtl_results_w_geno_means;
run;

data eqtl.results_w_onengut_credible_means;
   set eqtl_og_results_w_geno_means;
   onengut_index_snp_cat=tranwrd(onengut_index_snp_cat,', .' ,'');
run;


/* Export as Supplementary Tables */

 /**********************************************************************
 *   PRODUCT:   SAS
 *   VERSION:   9.4
 *   CREATOR:   External File Interface
 *   DATE:      18FEB16
 *   DESC:      Generated SAS Datastep Code
 *   TEMPLATE SOURCE:  (None Specified.)
 ***********************************************************************/
    data _null_;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    %let _EFIREC_ = 0;     /* clear export record count macro variable */
    file '/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/SuppTable1.csv'
delimiter=',' DSD DROPOVER lrecl=32767;
    if _n_ = 1 then        /* write column names or labels */
     do;
       put
          "Gene ID"
       ','
          "Multigene Flag (exonic regions)"
       ','
          "Autoimmune gene flag"
       ','
          "T1D candidate gene flag"
       ','
          "Feature type"
       ','
          "Feature ID"
       ','
          "SNP ID"
       ','
          "Significant splicing cis-eQTL flag"
       ','
          "CD4+ FDR P"
       ','
          "CD8+ FDR P"
       ','
          "CD19+ FDR P"
       ','
          "CD4+ HTZ/HMZ1 ratio"
       ','
          "CD4+ HMZ2/HMZ1 ratio"
       ','
          "CD4+ HMZ1 expression"
       ','
          "CD4+ HTZ expression"
       ','
          "CD4+ HMZ2 expression"
       ','
          "CD8+ HTZ/HMZ1 ratio"
       ','
          "CD8+ HMZ2/HMZ1 ratio"
       ','
          "CD8+ HMZ1 expression"
       ','
          "CD8+ HTZ expression"
       ','
          "CD8+ HMZ2 expression"
       ','
          "CD19+ HTZ/HMZ1 ratio"
       ','
          "CD19+ HMZ2/HMZ1 ratio"
       ','
          "CD19+ HMZ1 expression"
       ','
          "CD19+ HTZ expression"
       ','
          "CD19+ HMZ2 expression"
       ;
     end;
   set  EQTL.RESULTS_SUMMARY_TABLE_W_MEANS   end=EFIEOD;
       format gene_id $36. ;
       format flag_multigene best12. ;
       format flag_autoimmune_gene best12. ;
       format flag_diabetes_gene best12. ;
       format feature_type $4. ;
       format feature_id $2475. ;
       format snp_id $15. ;
       format flag_eqtl_sig best12. ;
       format CD4_FDR_P best12. ;
       format CD8_FDR_P best12. ;
       format CD19_FDR_P best12. ;
       format cd4_hmz1_htz_hmz1_ratio best12. ;
       format cd4_hmz2_hmz1_ratio best12. ;
       format cd4_hmz1_expression $30. ;
       format cd4_htz_expression $30. ;
       format cd4_hmz2_expression $30. ;
       format cd8_hmz1_htz_hmz1_ratio best12. ;
       format cd8_hmz2_hmz1_ratio best12. ;
       format cd8_hmz1_expression $30. ;
       format cd8_htz_expression $30. ;
       format cd8_hmz2_expression $30. ;
       format cd19_hmz1_htz_hmz1_ratio best12. ;
       format cd19_hmz2_hmz1_ratio best12. ;
       format cd19_hmz1_expression $30. ;
       format cd19_htz_expression $30. ;
       format cd19_hmz2_expression $30. ;
     do;
     EFIOUT + 1;
     put gene_id $ @;
     put flag_multigene @;
     put flag_autoimmune_gene @;
     put flag_diabetes_gene @;
     put feature_type $ @;
     put feature_id $ @;
     put snp_id $ @;
     put flag_eqtl_sig @;
     put CD4_FDR_P @;
     put CD8_FDR_P @;
     put CD19_FDR_P @;
     put cd4_hmz1_htz_hmz1_ratio @;
     put cd4_hmz2_hmz1_ratio @;
     put cd4_hmz1_expression $ @;
     put cd4_htz_expression $ @;
     put cd4_hmz2_expression $ @;
     put cd8_hmz1_htz_hmz1_ratio @;
     put cd8_hmz2_hmz1_ratio @;
     put cd8_hmz1_expression $ @;
     put cd8_htz_expression $ @;
     put cd8_hmz2_expression $ @;
     put cd19_hmz1_htz_hmz1_ratio @;
     put cd19_hmz2_hmz1_ratio @;
     put cd19_hmz1_expression $ @;
     put cd19_htz_expression $ @;
     put cd19_hmz2_expression $ ;
     ;
   end;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  if EFIEOD then call symputx('_EFIREC_',EFIOUT);
  run;


  /**********************************************************************
  *   PRODUCT:   SAS
  *   VERSION:   9.4
  *   CREATOR:   External File Interface
  *   DATE:      18FEB16
  *   DESC:      Generated SAS Datastep Code
  *   TEMPLATE SOURCE:  (None Specified.)
  ***********************************************************************/
     data _null_;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     %let _EFIREC_ = 0;     /* clear export record count macro variable */
     file '/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/SuppTable2.csv'
 delimiter=',' DSD DROPOVER lrecl=32767;
     if _n_ = 1 then        /* write column names or labels */
      do;
        put�
           "SNP ID"
        ','
           "Index SNPs"
        ','
           "Gene ID"
        ','
           "T1D candidate gene flag"
        ','
           "Feature ID"
        ','
           "Feature type"
        ','
           "Feature position"
        ','
           "Multigene flag (exonic regions)"
        ','
           "Significant splicing cis-eQTL flag"
        ','
           "CD4+ FDR P"
        ','
           "CD8+ FDR P"
        ','
           "CD19+ FDR P"
        ','
          "Credible SNPs (r2>0.8)"
       ','
          "Credible SNP IDs (r2>0.8)"
       ','
          "Credible SNPs (0.8>r2>0.5)"
       ','
          "Credible SNP IDs (0.8>r2>0.5)"
       ','
          "Credible SNPs (r2<0.5)"
       ','
          "Credible SNP IDs (r2<0.5)"
       ','
          "CD4+ HTZ/HMZ1 ratio"
       ','
          "CD4+ HMZ2/HMZ1 ratio"
       ','
          "CD4+ HMZ1 expression"
       ','
          "CD4+ HTZ expression"
       ','
          "CD4+ HMZ2 expression"
       ','
          "CD8+ HTZ/HMZ1 ratio"
       ','
          "CD8+ HMZ2/HMZ1 ratio"
       ','
          "CD8+ HMZ1 expression"
       ','
          "CD8+ HTZ expression"
       ','
          "CD8+ HMZ2 expression"
       ','
          "CD19+ HTZ/HMZ1 ratio"
       ','
          "CD19+ HMZ2/HMZ1 ratio"
       ','
          "CD19+ HMZ1 expression"
       ','
          "CD19+ HTZ expression"
       ','
          "CD19+ HMZ2 expression"
       ;
     end;
  set  EQTL.RESULTS_W_ONENGUT_CREDIBLE_MEANS   end=EFIEOD;
      format snp_id $15. ;
      format onengut_index_snp_cat $120. ;
      format gene_id $36. ;
      format flag_diabetes_gene best12. ;
      format feature_id $2475. ;
      format feature_type $4. ;
      format feature_pos $50. ;
      format flag_multigene best12. ;
      format flag_eqtl_sig best12. ;
      format CD4_FDR_P best12. ;
      format CD8_FDR_P best12. ;
      format CD19_FDR_P best12. ;
      format num_credible_snps_80perc best12. ;
      format credible_snp_r2_80 $600. ;
      format num_credible_snps_50perc best12. ;
      format credible_snp_r2_50 $600. ;
      format num_credible_snps_lowperc best12. ;
      format credible_snp_r2_low $600. ;
      format cd4_hmz1_htz_hmz1_ratio best12. ;
      format cd4_hmz2_hmz1_ratio best12. ;
      format cd4_hmz1_expression $30. ;
      format cd4_htz_expression $30. ;
      format cd4_hmz2_expression $30. ;
      format cd8_hmz1_htz_hmz1_ratio best12. ;
      format cd8_hmz2_hmz1_ratio best12. ;
      format cd8_hmz1_expression $30. ;
      format cd8_htz_expression $30. ;
      format cd8_hmz2_expression $30. ;
      format cd19_hmz1_htz_hmz1_ratio best12. ;
      format cd19_hmz2_hmz1_ratio best12. ;
      format cd19_hmz1_expression $30. ;
      format cd19_htz_expression $30. ;
      format cd19_hmz2_expression $30. ;
     do;
       EFIOUT + 1;
       put snp_id $ @;
       put onengut_index_snp_cat $ @;
       put gene_id $ @;
       put flag_diabetes_gene @;
       put feature_id $ @;
       put feature_type $ @;
       put feature_pos $ @;
       put flag_multigene @;
       put flag_eqtl_sig @;
       put CD4_FDR_P @;
       put CD8_FDR_P @;
       put CD19_FDR_P @;
       put num_credible_snps_80perc @;
       put credible_snp_r2_80 $ @;
       put num_credible_snps_50perc @;
       put credible_snp_r2_50 $ @;
       put num_credible_snps_lowperc @;
       put credible_snp_r2_low $ @;
       put cd4_hmz1_htz_hmz1_ratio @;
       put cd4_hmz2_hmz1_ratio @;
       put cd4_hmz1_expression $ @;
       put cd4_htz_expression $ @;
       put cd4_hmz2_expression $ @;
       put cd8_hmz1_htz_hmz1_ratio @;
       put cd8_hmz2_hmz1_ratio @;
       put cd8_hmz1_expression $ @;
       put cd8_htz_expression $ @;
       put cd8_hmz2_expression $ @;
       put cd19_hmz1_htz_hmz1_ratio @;
       put cd19_hmz2_hmz1_ratio @;
       put cd19_hmz1_expression $ @;
       put cd19_htz_expression $ @;
       put cd19_hmz2_expression $ ;
       ;
     end;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    if EFIEOD then call symputx('_EFIREC_',EFIOUT);
    run;


proc export data=eqtl.results_w_onengut_credible_means outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/SuppTable2.csv' dbms=csv replace;
run;

