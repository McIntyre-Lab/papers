
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*

for the 4568 transcripts detected in only elevated ozone in only 1 genotype, are the tpm values low??

input:
    (1) work.ele_all_cnt                        
        IDs transcripts detected on only 1 of the 5 genotypes
            from DD+transcript_overlap_across_geno_amb_only_ele_only.sas

    (2) pacbio.rsem_exp_matrix_tpm0_w_flags     
        contains TPM values
            from create_exp_matrix_4_supll_and_count_02amm.sas


mean TPM per library ranges for transcripts detected in only 1 genotype:  from 2.4 to 10.6

*/   





data ele_ones ;
set ele_all_cnt ;
where sum_across_geno = 1 ;
keep transcriptID ;
run ;

proc sort data = ele_ones ;
by transcriptID ;
proc sort data = pacbio.rsem_exp_matrix_tpm0_w_flags ;
by transcriptID ;
run ;

data ele2_ones ;
merge ele_ones (in=in1) pacbio.rsem_exp_matrix_tpm0_w_flags (in=in2) ;
by transcriptID ;
if in1 ;
run;

proc sort data = ele2_ones ;
by transcriptID ;
run;

%macro flips (geno) ;

proc transpose data = ele2_ones out = flip_&geno. ;
by transcriptID ;
var &geno._: ;
run ;

data flip2_&geno. ;
set flip_&geno. ;
rename _name_ = sample ;
rename col1 = TPM ;
run ;

proc sort data = flip2_&geno. ;
by sample ;

title "&geno. ";
proc means data = flip2_&geno. ;
var tpm ;
by sample ;
output out = ele_only_&geno. mean = mean ;
run;
%mend ;

%flips (B73) ;
%flips (C123) ;
%flips (Hp301) ;
%flips (Mo17) ;
%flips (NC338) ;


data all_ele ;
format sample $16. ;
set ele_only_: ;
run ;

proc sort data = all_ele ;
by mean ;
run ;

/* mean TPM per library ranges from 2.4 to 10.6 */






