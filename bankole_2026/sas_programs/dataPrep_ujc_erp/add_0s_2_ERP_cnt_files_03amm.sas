

/*
create design file for ERP_plus containing all combos 
*/

libname lmm_rmg "!MCLAB/lmm_rmg_dros_data/sasdata";
libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata";
libname lmm_rlr "!MCLAB/lmm_rlr_head_data/sasdata";
libname sex "!MCLAB/sex_specific_splicing/sasdata";



%macro flip (species, genome) ;

data &species._erp1 ;
set  sex.&species._erp_cnts_sumTR_stack;
ERP = scan(erp_plus, 1, "|");
keep geneID ERP_plus  ;
run ;

proc sort data =  &species._erp1 nodups;
by _all_ ;
run;


data &species._erp3 ;
set &species._erp1 ;
length sample $12. samples1 $12. samples2 $12. samples3 $12. samples4 $12. samples5 $12. samples6 $12. samples7 $12. samples8 $12. samples9 $12. samples10 $12. samples11 $12. samples12 $12. ;
array samples[12] $ ("&species._F_rep1", "&species._F_rep2", "&species._F_rep3", "&species._F_rep4", "&species._F_rep5", "&species._F_rep6", "&species._M_rep1", "&species._M_rep2", "&species._M_rep3", "&species._M_rep4", "&species._M_rep5", "&species._M_rep6") ;
do i = 1 to 12;
    sample = samples[i];
    output ;
end ;
drop samples: i ;
run;


proc sort data =  &species._erp3  ;
by sample geneID  ERP_plus  ;
run;

proc sort data = sex.&species._erp_cnts_sumTR_stack ;
by sample geneID  ERP_plus  ;
run ;

data &species._erp_w0 ;
merge sex.&species._erp_cnts_sumTR_stack &species._erp4  ;
by sample geneID  ERP_plus  ;
run;

data &species._erp_cnts_stack_with_0 ;
set  &species._erp_w0 ;
if ERP_sumReadNum = . then ERP_sumReadNum = 0;
else ERP_sumReadNum = ERP_sumReadNum ;
run ;

%mend ;

%flip (dmel, dmel6) ;
%flip (dsim, dsim2) ;
%flip (dser, dser1) ;


/* do for yak and san */
%macro flip2 (species, genome) ;

data &species._erp1 ;
set  sex.&species._erp_cnts_sumTR_stack;
ERP = scan(erp_plus, 1, "|");
keep geneID ERP  ;
run ;

proc sort data =  &species._erp1 nodups;
by _all_ ;
run;

data &species._erp2 ;
set &species._erp1 ;
    do flagDataOnlyExon = 0 to 1;
        do flagIR = 0 to 1;
            output ;
        end ;
    end ;
run;

data &species._erp3 ;
set &species._erp2 ;
length sample $12. samples1 $12. samples2 $12. ;
array samples[2] $ ("&species._F_rep1", "&species._M_rep1") ;
do i = 1 to 2;
    sample = samples[i];
    output ;
end ;
drop samples: i ;
run;

data &species._erp4 ;
set  &species._erp3 ;
ERP_plus = compress(ERP||'|'||flagDataOnlyExon||'|'||flagIR) ;
keep sample geneID ERP_plus ;
run ;

proc sort data =  &species._erp4  ;
by sample geneID  ERP_plus  ;
run;

proc sort data = sex.&species._erp_cnts_sumTR_stack ;
by sample geneID ERP_plus;
run ;

data &species._erp_w0 ;
merge sex.&species._erp_cnts_sumTR_stack &species._erp4  ;
by sample geneID ERP_plus;
run;

data &species._erp_cnts_stack_with_0 ;
set  &species._erp_w0 ;
if ERP_sumReadNum = . then ERP_sumReadNum = 0;
else ERP_sumReadNum = ERP_sumReadNum ;
run ;

%mend ;

%flip2 (dyak, dyak2) ;
%flip2 (dsan, dsan1) ;


%macro perms (species) ;

data sex.&species._erp_cnts_sumtr_stack_w0s ;
retain sample geneID ERP_plus ;
set &species._erp_cnts_stack_with_0 ;    
run ;
%mend ;

%perms (dmel) ;
%perms (dsim) ;
%perms (dser) ;
%perms (dyak) ;
%perms (dsan) ;


