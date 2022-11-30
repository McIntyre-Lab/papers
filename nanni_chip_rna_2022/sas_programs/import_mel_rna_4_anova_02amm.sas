
libname dros 'C:\a1stuff\dros';

libname dros "!MCLAB/Dros_PB_ChIP/sasdata/RNAseq";

libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';


/* prep for anova / ttests
    only keep apn_uq_ff data 
    drop offs
    merge in design
*/


/* keep only variables that end in uq_ff and featureID and on_off flags*/
%let suffix=_uq_ff;

proc sql noprint;
    select name into :keepList separated by ' '
    from dictionary.columns
    where upcase(libname)=upcase('chiprna')
    and upcase(memname)=upcase('mel_chip_rna_frag_flags_anno')
    and upcase(substr(name,length(name)-(length("&suffix")-1),length("&suffix")))=upcase("&suffix") ;
quit;

%put &keepList ;

data wide ;
set chiprna.mel_chip_rna_frag_flags_anno ( keep = &vars featureID flag_mel_m_on0_apn flag_mel_f_on0_apn) ;
run;

proc contents data = wide ; run;

/* drop offs */
data mel_flag_analyze;
set  wide ;
where flag_mel_m_on0_apn=1 or flag_mel_f_on0_apn=1;
drop 
flag_mel_m_on0_apn
flag_mel_f_on0_apn
;
run;

proc sort data = mel_flag_analyze ;
by featureID ;
run;

/*stack*/
proc transpose data=mel_flag_analyze out=mel_stack_4_anova;
by featureid;
run;
/* 3,496,152 or 145,673*24 --> has ALL fragments */


data mel_stack_4_anova2;
set mel_stack_4_anova;
rename _NAME_ = sampleID
col1=apn_uq_ff;
run;

/*LMM note - AMM will follow up:
need to recheck the ff  I don't like the discontinuity in how we did it*/

proc sort data=mel_stack_4_anova2;
by sampleid;
run;

proc import datafile = "!MCLAB/svn_lmm_dros_head_data/ethanol_srna/design_files/mel_rna_df_v2.csv"
out = mel_rna_df_v2 
dbms = csv replace ;
run;

proc sort data=mel_rna_df_v2;
by sampleid;
run;

/*note only saving perm here for convenience final version should be temp*/
data dTemp.mel_ready_4_anova oops ; 
merge mel_stack_4_anova2 (in=in1) mel_rna_df_v2 (in=in2);
by sampleid;
if in1 and in2 then output dTemp.mel_ready_4_anova ; 
else output oops ;
run;

