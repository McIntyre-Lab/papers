
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

libname anno "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/*      

create table containing flag for differentially detected genes where diff detected is on in amb OR on in ozone only
by genotype
       
    input
        pacbio.sub_geno_trt_&type._onCall_tpm0

    output
        pacbio.Diff_detect_amb_ele_gene 


**** checking what went to tappas: -- detected in both  

*/

title "" ;

%macro charts (genotype, type, ID) ;

data &genotype._&type. ;
format &ID $44. ;
set pacbio.sub_geno_trt_&type._onCall_tpm0;
keep flag_&genotype._: &ID ;
run;

data &genotype._&type._2 ;
set &genotype._&type. ;
length category $20. ;
length genotype $5. ;
if flag_&genotype._amb_on0 = 0 and flag_&genotype._ele_on0 = 0  then category = "Not Detected";
else if flag_&genotype._amb_on0 = 0 and flag_&genotype._ele_on0 = 1  then category = "Elevated Only";
else if flag_&genotype._amb_on0 = 1 and flag_&genotype._ele_on0 = 0  then category = "Ambient Only";
else if flag_&genotype._amb_on0 = 1 and flag_&genotype._ele_on0 = 1  then category = "Both Conditions";
else category = "oops" ;
drop flag_: ;
genotype = "&genotype" ;
run ;

data &genotype._&type._3 ;
set  &genotype._&type._2 ;
rename category = &genotype ;
if category = "Both Conditions" then flag_DD_&genotype = 0 ;
else if category = "Not Detected" then flag_DD_&genotype = 0 ;
else if category = "Elevated Only" then flag_DD_&genotype = 1 ;
else if category = "Ambient Only" then flag_DD_&genotype = 1 ;
else flag_DD_&genotype = 2 ;
drop genotype ;
run ;

%mend ;

%charts (B73, gene, geneID);
%charts (C123, gene, geneID);
%charts (Hp301, gene, geneID);
%charts (Mo17, gene, geneID);
%charts (NC338, gene, geneID);

data DD_all ;
merge B73_gene_3 C123_gene_3 Hp301_gene_3 Mo17_gene_3 NC338_gene_3 ;
by geneID ;
run;

proc freq data = DD_all ;
tables flag_DD_: ;
run ;

data pacbio.Diff_detect_amb_ele_gene ;
retain geneID DD_geno_count ;
set DD_all ;
DD_geno_count = sum(flag_DD_B73, flag_DD_C123, flag_DD_Hp301, flag_DD_Mo17, flag_DD_NC338) ; 
run ;





