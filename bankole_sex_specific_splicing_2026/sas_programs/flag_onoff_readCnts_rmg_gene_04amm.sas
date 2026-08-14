

libname lmm_axk "!MCLAB/lmm_axk_head_data/sasdata";
libname lmm_rmg "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sasdata" ;
libname sex "!MCLAB/sex_specific_splicing/sasdata";


/*

gene onoff

flag on/off using readNums of 0, 5, 10 and 25:  
    if off in F and off in M then   flag_both_sex_off_readNum&num
    if on in F and on on M then     flag_both_sex_on_readNum&num
    if off in F and on in M then    flag_Monly_on_readNum&num
    if on in F and off in M then    flag_Fonly_on_readNum&num

*** decided - go wth readCnt > 25 for GENE

for norm:
   drop if flag_both_sex_off_readNum10 = 1 

*/


data design_r ;
set lmm_rmg.design_rmg_mel_sim; 
drop sampleID ref techrep;
run;

proc sort data = design_r nodups ;
by _all_ ;
run;

data design_k ;
set lmm_axk.design_axk_ser; 
drop sampleID ref techrep;
run;

proc sort data = design_k nodups ;
by _all_ ;
run;




%macro import (species, ref, dsgn) ;

data &species.2&ref._gene2 ;
set sex.&species.2&ref._gene_cnts_from_jxnHash ;

if numTranscripts > 0 then readCnt_on0 = 1; else readCnt_on0 = 0;
if numTranscripts >= 5 then readCnt_on5 = 1; else readCnt_on5 = 0;
if numTranscripts >= 10 then readCnt_on10 = 1; else readCnt_on10 = 0;
if numTranscripts >= 25 then readCnt_on25 = 1; else readCnt_on25 = 0;
run;

proc sort data = &dsgn;
by sample ;
run ;

data &species.2&ref._gene_df ;
merge &species.2&ref._gene2 (in=in1) &dsgn (in=in2);
by sample ;
if in1 ;
run;

data &species.2&ref._gene_df2 ;
set &species.2&ref._gene_df ;
where rep = "rep4" or rep = "rep5" or rep = "rep6";
run ; 

%mend ;

%import (mel, dmel6, design_r) ;
%import (sim, dsimW, design_r) ;
%import (sim, dsim2, design_r) ;
%import (ser, dser1, design_k) ;



%macro which (species, ref) ;
%macro oncalls (num) ;

proc sort data=&species.2&ref._gene_df2 ;
  by sex;
  run;

/* Going to use readNum 0, 5 and 10 for flagging on/off */
proc sort data= &species.2&ref._gene_df2 ;
by sex geneID ;
run;


/* gene is on if expressed in in 50% of the reps */
proc means data = &species.2&ref._gene_df2  noprint ;
    by sex geneID ;
    var readCnt_on&num;
    output out = &species.2&ref._sex_on&num mean=sex_percent_on ;
    run ;

proc sort data = &species.2&ref._sex_on&num ;
    by geneID ;
    run;

proc transpose data = &species.2&ref._sex_on&num out = &species.2&ref._sex_on&num._sbys ;
    by geneID ;
    id sex ;
    var sex_percent_on ;
    run;

data sex.onCalls_gene_&species.2&ref._RC&num..;              * if 50% of reps then is expressed ;
    set &species.2&ref._sex_on&num._sbys ;                    * using sex to determine ;
    by geneID ;

    if  f > 0.5 then flag_F_on = 1 ; else flag_F_on = 0 ;
    if  m > 0.5 then flag_M_on = 1 ; else flag_M_on = 0 ;

    /* if off in F and off in M then flag all off */
    if flag_F_on = 0 and flag_M_on = 0 then flag_both_sex_off_readCnt&num = 1 ;
    else flag_both_sex_off_readCnt&num = 0;

    /* if on in F and on in M then flag all on */
    if flag_F_on = 1 and flag_M_on = 1 then flag_both_sex_on_readCnt&num = 1 ;
    else flag_both_sex_on_readCnt&num = 0 ;

    /* if off in F and on in M then flag M only on */
    if flag_F_on = 0 and flag_M_on = 1 then flag_Monly_on_readCnt&num = 1;
    else flag_Monly_on_readCnt&num = 0 ;

    /* if on in F and off in M then flag F only on */
    if flag_F_on = 1 and flag_M_on = 0 then flag_Fonly_on_readCnt&num = 1 ;
    else flag_Fonly_on_readCnt&num = 0 ;
    keep geneID flag_: f m ;
    run ;
%mend ;

%onCalls (0) ;
%onCalls (5) ;
%onCalls (10) ;
%onCalls (25) ;
%mend ;
%which (mel, dmel6) ;
%which (sim, dsimW) ;
%which (sim, dsim2) ;



ods pdf file = "!MCLAB/sex_specific_splicing/expression_output/freqs_onCalls_gene_mel_sim_ser.pdf" ;

%macro which (species, ref) ;
%macro freqs (num) ;

title "GENE &species.2&ref MF on calls,  readCnt > &num";
proc freq data = sex.onCalls_gene_&species.2&ref._RC&num.. ;
tables flag_:  ;
run;

proc freq data =sex.onCalls_gene_&species.2&ref._RC&num.. ;
tables flag_F_on * flag_M_on  ;
run;

title "";
%mend ;

%freqs (25) ;
%freqs (0) ;
%freqs (5) ;
%freqs (10) ;
%mend ;
%which (mel, dmel6) ;
%which (sim, dsimW) ;
%which (sim, dsim2) ;

ods pdf close ;




