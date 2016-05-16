filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/cegs_ase_paper/sas_data';
libname cegs2 'Z:/cegs_ase_paper/sas_data';
libname svn 'Z:/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/sas_data';

 data WORK.ase_genes    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile 'Z:\cegs_ase_paper\output\ase_gene_list.csv' delimiter =
 ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
          informat symbol $8. ;
          informat gene_compare_enrich $8. ;
          format symbol $8. ;
          format gene_compare_enrich $8. ;
       input
                   symbol $
                   gene_compare_enrich $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
	run;

data flag_ase_genes;
set ase_genes;
flag_ase_gene=1;
run;

proc sort data=cegs.cis_calls;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data cis_calls_w_gene;
merge cegs.cis_calls (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

data pool_genes;
set cegs.pool_2015_gene_lists;
rename gene_name = symbol;
run;

data rita_genes;
set cegs.rita_results_w_flags;
rename gene_symbol=symbol;
run;

proc sort data=rita_genes;
by symbol;
proc sort data=pool_genes;
by symbol;
proc sort data=flag_ase_genes;
by symbol;
run;

data combine_lists;
merge rita_genes pool_genes flag_ase_genes (in=in1);
by symbol;
if in1;
run; *138830 obs;

data chk_1;
set combine_lists;
if mel_biased=1 and flag_ase_gene=1;
run;
*flags to consider: sim_biased, mel_biased, our_defense, MK_adaptive_any, flag_table4_african_gene,
flag_table4_europe_gene, flag_table7_gene;

proc freq data=combine_lists;
by symbol;
tables flag_ase_gene*mel_biased / chisq;
ods output chiSq=mel_bias_chi_sq_output;
run;

data mel_bias_chi_sq_output2;
set mel_bias_chi_sq_output;
if statistic ne "Chi-Square" then delete;
run;
