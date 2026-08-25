libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata" ;
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata" ;
libname sex "!MCLAB/sex_specific_splicing/sasdata" ;


/* 


normalize using gene norm FF

jxnHash var
 

dmel, dsim, dser 
reps 4-6 only


*/

/* mel, sim, ser*/

%macro ratio (species, var ) ;

data &species._ujc_cnts ;
set sex.&species._ujc_cnts_sumtr_stack_w0s ;
    sex = scan(sample, 2, '_') ;
    rep = scan(sample, 3, '_') ;
    species = scan(sample, 1, '_') ;
run;

proc sort data = &species._ujc_cnts ; 
by jxnHash ;
proc sort data = sex.&species._jxnhash2geneid ;
by jxnHash ;
run;

data &species._ujc_cnts2 ;
merge &species._ujc_cnts   sex.&species._jxnhash2geneid ;
by jxnHash ;
run;


proc sort data =  &species._ujc_cnts2  ;
by geneID jxnHash ;
proc sort data = sex.onCall_ujc_&species._on0 ;
by geneID jxnHash ;
run ;

data &species._ujc_cnts_df ;
merge  &species._ujc_cnts2 (in=in1) sex.onCall_ujc_&species._on0  (in=in2) ;
by geneID jxnHash ;
run;

data ujc_&species. ;
set &species._ujc_cnts_df;
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

proc sort data = ujc_&species.  ;
by sample ;
proc sort data = ff_&species. ;
by sample ;
run;

data ff_ujc_&species. ;
merge ff_&species. (in=in1) ujc_&species. (in=in2) ;
by sample ;
if in2  ;
run ;

data ff2_ujc_&species. ;
retain ujc sample log_uq_sumReadNum ;
set ff_ujc_&species. ;
log_uq_sumReadNum = (uq_ff * jxnHash_sumReadNum) ;
run ;

data sex.norm_ujc_&species. ;
set ff2_ujc_&species. ;
run ;

%mend ;

%ratio (dmel, 762.5) ; 

%ratio (dsim, 898 ) ;
%ratio (dser, 752) ;





