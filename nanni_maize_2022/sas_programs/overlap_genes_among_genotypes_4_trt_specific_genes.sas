
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

libname anno "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/*      Do the genes/isoforms detected in only ambient or only elevated for each genotype overlap?   --> differentially detected 
    
for each genotype 
    identified genes detected in only elevated or only ambient conditions
   
combine genotypes
    for ambient only genes - count number of genotypes detected in
    for elevated only genes - count number of genotypes detected in

add geneName based on geneID 
    
    input
        pacbio.sub_geno_trt_&type._onCall_tpm0
        anno.B73_V4_geneID_geneName

    output
        MCLAB/maize_ozone_FINAL/pacbio_paper/figs/Number_genes_transcripts_detected_in_oz_or_amb_only_shared_among_genotypes.pdf


**** checking what went to tappas: -- detected in both  

*/



%macro charts (genotype, type, ID) ;

data &genotype._&type. ;
length &ID $44. ;
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

data &genotype._ele_&type. &genotype._amb_&type. ;
set  &genotype._&type._2 ;
rename category = &genotype ;
if find(category, "Elevated") ge 1 then output &genotype._ele_&type. ;
if find(category, "Ambient") ge 1  then output &genotype._amb_&type. ;;
drop genotype ;
run ;

proc sort data = &genotype._ele_&type. ;
by &ID ;
proc sort data = &genotype._amb_&type. ;
by &ID ;
run;

%mend ;

%charts (B73, gene, geneID);
%charts (C123, gene, geneID);
%charts (Hp301, gene, geneID);
%charts (Mo17, gene, geneID);
%charts (NC338, gene, geneID);

%charts (B73, isoform, transcriptID);
%charts (C123, isoform, transcriptID);
%charts (Hp301, isoform, transcriptID);
%charts (Mo17, isoform, transcriptID);
%charts (NC338, isoform, transcriptID);


ods pdf file = "!MCLAB/maize_ozone_FINAL/pacbio_paper/figs/Number_genes_transcripts_detected_in_oz_or_amb_only_shared_among_genotypes.pdf" ;

%macro nextBit (type, ID) ;

data all_ele_&type;
length &ID $ 44.;
merge B73_ele_&type C123_ele_&type Hp301_ele_&type Mo17_ele_&type NC338_ele_&type ;
by &ID ;
run ;

data all_ele_&type.2 ;
set all_ele_&type. ;
where B73  ne "" or Mo17  ne "" or C123  ne "" or Hp301  ne "" or NC338  ne "" ; 
count = 5 - cmiss(b73,  c123, hp301, mo17, nc338) ; 
run;   

proc sort data = all_ele_&type.2 ;
by &ID ;
run ;

title "Num &type. detected in elevated-only, shared among genotypes";
proc freq data = all_ele_&type.2 ;
tables count ;
run;

data all_amb_&type. ;
length &ID $ 44.;
merge B73_amb_&type. C123_amb_&type. Hp301_amb_&type. Mo17_amb_&type. NC338_amb_&type. ;
by &ID ;
run ;

data all_amb_&type.2 ;
set all_amb_&type. ;
where B73  ne "" or Mo17  ne "" or C123  ne "" or Hp301  ne "" or NC338  ne "" ; 
count = 5 - cmiss(b73,  c123, hp301, mo17, nc338) ; 
run;   

proc sort data = all_amb_&type.2 ;
by &ID ;
run ;


title "Num &type. detected in ambient-only, shared among genotypes";
proc freq data = all_amb_&type.2 ;
tables count ;
run;

%mend ;

%nextBit (gene, geneID);
%nextBit (isoform, transcriptID);

ods pdf close ;



/* merge in geneName  -- geneID */


data wanno_4_genes;
set pacbio.fsm_ism_nic_nnc_wanno;
drop isoform ;
run ;

proc sort data = wanno_4_genes nodups ;
by _all_ ;
run;

proc sort data = wanno_4_genes;
by geneID ;
run ;

%macro GA (trt) ;

proc sort data = all_&trt._gene2 ;
by geneID ;
run;

data all_&trt._gene3 ;
merge  all_&trt._gene2 (in=in1) wanno_4_genes (in=in2) ;
by geneID ;
if in1 ;
run;

data all_&trt._gene_ready ;
retain geneID ZMgn geneName count;
set  all_&trt._gene3 ;
rename count = number_genotypes_detected ;
run ;

proc export data = all_&trt._gene_ready
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/supplementary_data/overlap_genes_&trt._only_among_genotypes.txt"
dbms =tab replace ;
run;

%mend ;

%GA (ele) ;
%GA (amb) ;


data wanno_4_isoforms;
set pacbio.fsm_ism_nic_nnc_wanno;
rename isoform = transcriptID ;
run ;

proc sort data = wanno_4_isoforms nodups ;
by _all_ ;
run;

proc sort data = wanno_4_isoforms;
by transcriptID ;
run ;


%macro TA (trt) ;

proc sort data = all_&trt._isoform2;
by transcriptID ;
run;

data all_&trt._isoform3 ;
merge  all_&trt._isoform2 (in=in1) wanno_4_isoforms (in=in2) ;
by transcriptID ;
if in1 ;
run;

data all_&trt._iso_ready ;
retain transcriptID geneID ZMgn geneName count;
set  all_&trt._isoform3 ;
rename count = number_genotypes_detected ;
run ;

proc export data = all_&trt._iso_ready
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/supplementary_data/overlap_isoforms_&trt._only_among_genotypes.txt"
dbms =tab replace ;
run;

%mend ;

%TA (ele) ;
%TA (amb) ;

/* make perm */

data pacbio.diff_detect_ele_isoform ;
set all_ele_iso_ready ;
run ;

data pacbio.diff_detect_amb_isoform ;
set all_amb_iso_ready ;
run ;
data pacbio.diff_detect_ele_gene ;
set all_ele_gene_ready ;
run ;
data pacbio.diff_detect_amb_gene ;
set all_amb_gene_ready ;
run ;


proc sort data = all_ele_iso_ready ;
by descending number_genotypes_detected ;
run;

