libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata" ;
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata" ;
libname sex "!MCLAB/sex_specific_splicing/sasdata" ;


/* 


normalize using gene norm FF

ERP_plus var
 

dmel, dsim, dser 
reps 4-6 only


*/




/* mel, sim, ser*/

%macro ratio (species, var ) ;

data &species._erp_cnts ;
set sex.&species._erp_cnts_sumtr_stack_w0s ;
    sex = scan(sample, 2, '_') ;
    rep = scan(sample, 3, '_') ;
    species = scan(sample, 1, '_') ;
run;


proc sort data =  &species._erp_cnts  ;
by geneID ERP_plus ;
proc sort data = sex.onCall_ERP_&species._on0 ;
by geneID ERP_plus ;
run ;

data &species._ERP_cnts_df ;
merge  &species._erp_cnts (in=in1) sex.onCall_ERP_&species._on0  (in=in2) ;
by geneID ERP_plus ;
run;

data ERP_&species. ;
set &species._ERP_cnts_df;
where rep = "rep4" or rep = "rep5" or rep = "rep6";
if flag_both_sex_off_readCnt0 = 0 ;
run ; 


data ff_&species. ;
retain sample ;
set sex.&species._q3_med_totals ;
        uq_ff = &var./q3;
sample = compress("&species._"||sex||'_'||rep) ;
keep uq_ff sample ;
run;
        
/* merge mel_FF to mel_ts by sample 
   log_ug_readCnt = (uq_ff * FSM_ISM_count)

 */

proc sort data = ERP_&species.  ;
by sample ;
proc sort data = ff_&species. ;
by sample ;
run;

data ff_ERP_&species. ;
merge ff_&species. (in=in1) ERP_&species. (in=in2) ;
by sample ;
if in2  ;
run ;

data ff2_ERP_&species. ;
retain ERP sample log_uq_sumReadNum ;
set ff_ERP_&species. ;
log_uq_sumReadNum = (uq_ff * ERP_sumReadNum) ;
run ;

data sex.norm_ERP_&species. ;
set ff2_ERP_&species. ;
run ;

%mend ;

%ratio (dmel, 762.5) ; 

%ratio (dsim, 898 ) ;
%ratio (dser, 752) ;





