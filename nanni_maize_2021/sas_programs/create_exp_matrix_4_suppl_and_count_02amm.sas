libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";


/*

create expression matrix of data going into tappas including onCall flags

    input:  pacbio.sub_geno_trt_isoform_oncall_tpm0 
            pacbio.sbys_&geno._4_tappas
            
    
    output: pacbio.rsem_exp_matrix_tpm0_w_flags     
            MCLAB/maize_ozone_FINAL/pacbio_paper/penultimate_version/rsem_exp_matrix_tpm0_w_flags.tsv

            flag_&genom._into_tappas = 1 if in pacbio.sbys_&geno._4_tappa   
        
            ***NOTE:  flag_&geno._into_tappas = flag_&geno._trt_trnscpt_on0 in output file  (on in both conditions!!)

    # transcripts and genes into on-off calls (after rsem) 
        33,267 transcripts and 12,702 genes

    # transcripts and genes into tappas (i.e. "on" in both conditions)
        30,302 transcripts on in at least 1 genotype
        12,680 genes on in at least 1 genotype 



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

data ready2_&geno. ;
set ready_&geno. ;
*keep transcriptID tpm_noZero sample ;
*rename tpm_noZero = tpm ;
run;

/* create sbys */
proc transpose data = ready2_&geno out = flip_&geno._sbys ;
by transcriptID ;
id sample ;
var tpm ;
run;

data sbys_&geno._4_EM;
set  flip_&geno._sbys;   
drop _name_ ;
run;

proc sort data = sbys_&geno._4_EM; 
by transcriptID ;
run;

%mend ;

%separate_geno (B73) ;
%separate_geno (Mo17) ;
%separate_geno (C123) ;
%separate_geno (NC338) ;
%separate_geno (Hp301) ;

data sbys_4_em ;
merge sbys_: ;
by transcriptID ;
run;


data oncalls;
set pacbio.sub_geno_trt_isoform_oncall_tpm0 ;
run ;

data exp_matrix ;
merge sbys_4_em oncalls ;
by transcriptID ;
run;

data exp_matrix_2 ;
set  exp_matrix;
if flag_B73_trt_trnscpt_on0 = 1 or flag_C123_trt_trnscpt_on0 = 1 or flag_Hp301_trt_trnscpt_on0 = 1 or flag_Mo17_trt_trnscpt_on0 = 1 or flag_NC338_trt_trnscpt_on0 = 1 
    then flag_into_tappas = 1 ;
else flag_into_tappas = 0;
run ;

proc freq data = exp_matrix_2 ;
tables flag_into_tappas ;
run;   /* 30,302 transcripts into tappas */

data pacbio.rsem_exp_matrix_tpm0_w_flags ;
set exp_matrix_2 ;
run;


proc export data = pacbio.rsem_exp_matrix_tpm0_w_flags  
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/penultimate_version/rsem_exp_matrix_tpm0_w_flags.tsv"
dbms = tab replace ;
run ;




