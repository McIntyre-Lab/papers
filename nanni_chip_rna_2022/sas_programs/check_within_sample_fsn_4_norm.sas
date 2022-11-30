libname dros "!MCLAB/Dros_PB_ChIP/sasdata/RNAseq";


/*
fusions

LMM and AMM looked at mel vs sim and noEtoh vs etoh 
    
    *** median Q3 very similar, therefore selecting 63 for calculating uq_ff (in normalize_log_uq3_02amm.sas)
*/



%macro norm (species, trt) ;

/* merge in design file */
proc sort data = dros.df_pbrna_4_analysis ;
by sampleID ;
proc sort data = dros.cnts_&species._fusions_apn_stack ;
by sampleID ;
run;

data &species._fusions ouch_&species;
merge dros.cnts_&species._fusions_apn_stack (in=in1) dros.df_pbrna_4_analysis (in=in2) ;
by sampleID ;
if in1 and in2 then output &species._fusions ;
else output ouch_&species. ;
run;


/* Create clean dataset - drop fusions with little expression (keep flag_fusion_F_on = 1 and flag_fusion_M_on = 1) */
proc sort data=&species._fusions;
by featureID ;
proc sort data = dros.onCalls_fsn_&species._etoh_gt_apn0 ;
by featureID ;

proc sort data = dros.onCalls_fsn_&species._noEtoh_gt_apn0 ;
by featureID ;
run;

data &species._clean1;
merge &species._fusions (in=in1) dros.onCalls_fsn_&species._etoh_gt_apn0  ;
by featureID ;
if flag_fusion_F_on = 1 and flag_fusion_M_on = 1;
drop flag_: ;
run;  

data &species._clean ;
merge &species._clean1 (in=in1) dros.onCalls_fsn_&species._noEtoh_gt_apn0 ;
by featureID ;
if flag_fusion_F_on = 1 and flag_fusion_M_on = 1;
drop flag_: ;
run;  

/* Calculate With-in sample basic statistics */
proc means data=&species._clean noprint;
    class sampleID ;
    var reads_in_region;
    output out=&species._mapped sum=sum_mapped q1=q1 q3=q3 median=median mean=mean;
    run;

data &species._mapped_reads  ;
set &species._mapped ;
where _type_ = 1 ;
drop _type_ _freq_ ;
run ;
%mend ;

%norm (mel) ;
%norm (sim ) ;

data dros.check_fsn_4_norm ;
set mel_mapped_reads sim_mapped_reads ;
run ;

proc export data = dros.check_fsn_4_norm 
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/normalization/check_fsn_4_norm.tsv"
dbms = tab replace ;
run ;


ods pdf file = "!MCLAB/Dros_PB_ChIP/RNAseq/normalization/check_within_sample_fusions_4_norm.pdf" ;

proc plot data =  dros.check_fsn_4_norm  ;
plot sum_mapped * q3 ;
run;

proc sort data =  dros.check_fsn_4_norm  ;
by q3 ;
run;   /* sim_12_m_noEtoh_rep2 has low q3!! */

proc rank data =  dros.check_fsn_4_norm out = ranked ;
var  sum_mapped q3 ;
run;

title "all" ;
proc means data = dros.check_fsn_4_norm  q1 q3 mean median min max ;
run ;

data all_mel ;
set dros.check_fsn_4_norm ;
if find(sampleID, "mel") ge 1 ;
run;

data all_sim ;
set dros.check_fsn_4_norm ;
if find(sampleID, "sim") ge 1 ;
run;
title "all mel" ;
proc means data = all_mel  q1 q3 mean median min max ;
run ;

title "all sim" ;
proc means data = all_sim  q1 q3 mean median min max ;
run ;


data noEtoh ;
set dros.check_fsn_4_norm ;
if find(sampleID, "noEtoh") ge 1 ;
run;

title "noEtoh" ;
proc means data = noEtoh q1 q3 mean median min max;
run ;

data Etoh ;
set dros.check_fsn_4_norm ;
if find(sampleID, "_etoh") ge 1 ;
run;

title "Etoh" ;
proc means data = noEtoh q1 q3 mean median min max;
run ;

title "";

ods pdf close ;














