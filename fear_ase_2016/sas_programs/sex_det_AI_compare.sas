filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/cegs_ase_paper/sas_data';
libname svn 'Z:/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/sas_data';


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

 data WORK.sex_det    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile 'Z:\dros_human_sem_links\gene_lists\sex_det.csv' delimiter = ',' MISSOVER DSD
 lrecl=32767 firstobs=2 ;
         informat symbol $6. ;
         informat FBgn $11. ;
         format symbol $6. ;
         format FBgn $11. ;
      input
                  symbol $
                  FBgn $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

proc sort data=sex_det;
by symbol;
proc sort data=cis_calls_w_gene;
by symbol;
run;
