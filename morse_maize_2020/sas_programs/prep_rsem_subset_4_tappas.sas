


/*
prep for tappas
    data for each genotype
        data should have missing values, not zero's
        drop if flag_&geno._trt_trnscpt_on0 = 0
        
note - bad libraries:   B73_P1_C7_Ele
                        B73_P4_C1_Amb
                        B73_p4_C6_Amb                    
*/

data DF ;
set pacbio.design_file_no_failed_libs ;
new1 = tranwrd(sample, "Hp30_", "Hp301_") ;
new2 = tranwrd(new1, "NC33_", "NC338_") ;
geno1 = tranwrd(genotype, "NC333", "NC338") ;
drop sample genotype new1 ;
run ;

data design ;
set df ;
rename new2 = sample ;
rename geno1 = genotype ;
run ;

/* merge in design file */
proc sort data = pacbio.rsem_subset_isoforms_tpm_stk ;
by sample ;
proc sort data = design ;
by sample ;
run;

data datain oops ;
merge pacbio.rsem_subset_isoforms_tpm_stk (in=in1) design (in=in2) ; 
by sample;
if in1 and in2 then output datain ;
else output oops ;
run;

/* make 0's missing */
data datain2 ;
set datain ;
if tpm = 0 then tpm_noZero = .;
    else tpm_noZero = tpm ;
run ;

data datain2 ;
retain transcriptID sample tpm tpm_noZero ;
set datain2 ;
run;

/* merge in on/off calls */
proc sort data =pacbio.sub_geno_trt_isoform_oncall_tpm0 ;
by transcriptID ;
proc sort data = datain2 ;
by transcriptID ;
run ;

data with_flags oops ;
merge datain2 (in=in1) pacbio.sub_geno_trt_isoform_oncall_tpm0 (in=in2);
by transcriptID ;
if in1 and in2 then output with_flags ;
else output oops ;
run;  /* 0 in oops */

data pacbio.exp_matrix_tpm0_w_onCalls_stk ;
set with_flags ;
run ;


/* split out genotypes */
%macro separate_geno (geno) ;

data split_&geno. ;
set datain2 ;
where genotype = "&geno" ;
run;

/* merge in on/off calls */
proc sort data = split_&geno ;
by transcriptID ;
proc sort data = pacbio.sub_geno_trt_isoform_oncall_tpm0 ;
by transcriptID ;
run ;

data ready_&geno. oops;
merge split_&geno (in=in1) pacbio.sub_geno_trt_isoform_oncall_tpm0 (in=in2);
by transcriptID ;
if in1 and in2 then output ready_&geno;
else output oops ;
run;

/* drop if flag_&geno._trt_trnscpt_on0 = 1 */
data flag_&geno. ;
set ready_&geno ;
where flag_&geno._trt_trnscpt_on0 = 1 ;
run ;

data tappas_ready_&geno. ;
set flag_&geno. ;
keep transcriptID tpm_noZero sample ;
rename tpm_noZero = tpm ;
run;
/* create sbys for tappas */
proc transpose data = tappas_ready_&geno out = tappas_&geno._sbys ;
by transcriptID ;
id sample ;
var tpm ;
run;

data pacbio.sbys_&geno._4_tappas ;
set  tappas_&geno._sbys;   
drop _name_ ;
run;
/*
proc export data = pacbio.sbys_&geno._4_tappas 
outfile = "!MCLAB//maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/sbys_&geno._4_tappas.tsv"
label
dbms = tab replace ;
run ;
*/
/* create df for each geno */
data df2_&geno. ;
set ready_&geno ;
keep sample ;
run;

proc sort data = df2_&geno. nodups ;
by _all_ ;
run ;

data pacbio.df_&geno._4_tappas ;
set df2_&geno. ;
label   sample = 'Sample'   
        condition = 'Condition'  ;
trt = compress(scan(sample, -1, '_')) ;
if trt = "Amb" then condition = "Ambient" ;
else if trt = "Ele" then condition = "Ozone" ;
else condition = "oops" ;
drop trt ;
run;
/*
proc export data = pacbio.df_&geno._4_tappas
outfile = "!MCLAB//maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/df_&geno._4_tappas.tsv"
label
dbms = tab replace ;
run ;
*/
%mend ;

%separate_geno (B73) ;
%separate_geno (Mo17) ;
%separate_geno (C123) ;
%separate_geno (NC338) ;
%separate_geno (Hp301) ;










