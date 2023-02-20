
libname dros 'C:\a1stuff\dros';

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';

libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";

proc contents data = chiprna.sim_chip_rna_frag_flags_anno ;run;

/* keep only variables that end in uq_ff and featureID and on_off flags*/
%let suffix=_uq_ff;

proc sql noprint;
    select name into :keepList separated by ' '
    from dictionary.columns
    where upcase(libname)=upcase('chiprna')
    and upcase(memname)=upcase('sim_chip_rna_frag_flags_anno')
    and upcase(substr(name,length(name)-(length("&suffix")-1),length("&suffix")))=upcase("&suffix") ;
quit;

%put &keepList ;

data sim_wide ;
set chiprna.sim_chip_rna_frag_flags_anno ( keep = &keepList featureID flag_sim_m_on0_apn flag_sim_f_on0_apn) ;
run;

proc contents data = sim_wide ; run;

/* drop offs and other extraneous variables */
data sim_flag_analyze;
set  sim_wide ;
where flag_sim_m_on0_apn=1 or flag_sim_f_on0_apn=1;
drop 
flag_: 
ratio_avg_f_m_apn_uq_ff
avg_: 
;
run;

proc sort data = sim_flag_analyze ;
by featureID ;
run;
proc contents data = sim_flag_analyze ; run ;

/*stack*/
proc transpose data=sim_flag_analyze out=sim_stack_4_anova;
by featureid;
run;
/*  23 * 135452 = 3115396 --> has ALL fragments */

proc freq data = sim_stack_4_anova ;
tables _name_ / out = cnt_sim_name ;
run;

data sim_stack_4_anova2;
set sim_stack_4_anova;
rename _NAME_ = sampleID
col1=apn_uq_ff;
run;

/*LMM note - AMM will follow up:
need to recheck the ff  I don't like the discontinuity in how we did it*/

proc sort data=sim_stack_4_anova2;
by sampleid;
run;

proc import datafile = "!MCLAB/svn_lmm_dros_head_data/ethanol_srna/design_files/sim_rna_df_v2.csv"
out = sim_rna_df_v2 
dbms = csv replace ;
run;

proc sort data=sim_rna_df_v2;
by sampleid;
run;

/*note only saving perm here for convenience final version should be temp*/
data dTemp.sim_ready_4_anova oops ; 
merge sim_stack_4_anova2 (in=in1) sim_rna_df_v2 (in=in2);
by sampleid;
if in1 and in2 then output dTemp.sim_ready_4_anova ; 
else output oops ;
run;

