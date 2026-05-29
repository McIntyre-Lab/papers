libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;


/*

flagAnalyze_frag_DE = 1

model DE by cell type
  by fragment 
  t1d sex age 
 
roll to gene 

note - normalized by cell type, kept frag if on in more than 1 group, kept frag if ge 10bp

*/



/*  prep data */
%macro prepping (cell) ;

/* add design */
data add_meta ;
set seq.illumina_RO1_meta_design;
keep sampleID enrollment_age population flag_female flag_case ;
run ;
proc sort data = add_meta ;
by sampleID ;
proc sort data = seq.norm_data_frag_&cell._stack  ;
by sampleID ;
run;

data almost_&cell. ;
merge seq.norm_data_frag_&cell._stack  (in=in1) add_meta (in=in2) ;
by sampleID  ;
if in1 ;
run;

data ready_DE_&cell. ;
set almost_&cell ;
keep featureID sampleID population flag_female flag_case log_uq_apn enrollment_age  ;
run ;

proc sort data = ready_DE_&cell ;
by featureID ;
proc sort data = seq.flag_analyze_frag_DE_&cell. ;
by featureID ;
run ;

data near_&cell. oops_&cell. ;
merge seq.flag_analyze_frag_DE_&cell. (in=in1) ready_DE_&cell  (in=in2) ;
by featureID   ;
if in2 then output near_&cell. ;
else output oops_&cell. ;
run;

data ready_DE_&cell. ;
set near_&cell. ;
where flag_analyze_frag_DE = 1 ;
run;

%mend ;

%prepping (CD4) ;
%prepping (CD8) ; 
%prepping (CD19) ; 



proc contents data = ready_DE_CD4 ; run;


%macro de2_model (cell) ;

proc glimmix data = ready_DE_&cell. ;
by featureID ;
class flag_female flag_case  ;
model log_uq_apn = enrollment_age flag_female flag_case / htype =1 ;
lsmeans flag_female flag_case ;
output out = seq.de_resid_frag_&cell._all resid=resid pred=pred student=student ;
ods output 
    tests1 = seq.de_t1_frag_&cell._all  
    lsmeans = seq.de_ls_frag_&cell._all ;
run;

%mend ;

%de2_model (CD4) ;
%de2_model (CD8) ;
%de2_model (CD19) ;

%macro de3_model (cell) ;

proc glimmix data = ready_DE_&cell. ;
by featureID ;
class flag_female flag_case  ;
model log_uq_apn = enrollment_age flag_female flag_case flag_female * flag_case / htype =1 ;
lsmeans flag_female * flag_case / slice = flag_female pdiff = all ;
output out = seq.de_resid_frag_&cell._FbyCase resid=resid pred=pred student=student ;
ods output 
    tests1 = seq.de_t1_frag_&cell._FbyCase 
    lsmeans = seq.de_ls_frag_&cell._FbyCase 
    diffs = seq.de_pd_frag_&cell._FbyCase  
    slices = seq.de_sl_frag_&cell._FbyCase  ;
run;

%mend ;

%de3_model (CD4) ;
%de3_model (CD8) ;
%de3_model (CD19) ;


%macro output (cell) ;

proc export data = seq.de_resid_frag_&cell._FbyCase
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex/de_resid_frag_&cell._FbyCase.csv"
dbms = csv replace ;
run;

proc datasets library = SEQ ;
modify de_t1_frag_&cell._FbyCase  ;
format probF pvalue32.30 ;
run;
proc export data = seq.de_t1_frag_&cell._FbyCase   
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex/de_t1_frag_&cell._FbyCase.csv"
dbms = csv replace ;
run;

proc datasets library = SEQ ;
modify de_ls_frag_&cell._FbyCase  ;
format probt pvalue32.30 ;
run;
proc export data = seq.de_ls_frag_&cell._FbyCase  
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex/de_ls_frag_&cell._FbyCase.csv"
dbms = csv replace ;
run;

proc datasets library = SEQ ;
modify seq.de_pd_frag_&cell._FbyCase  ;
format probt pvalue32.30 ;
run;
proc export data = seq.de_pd_frag_&cell._FbyCase  
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex/de_pd_frag_&cell._FbyCase.csv"
dbms = csv replace ;
run;

proc datasets library = SEQ ;
modify seq.de_sl_frag_&cell._FbyCase;
format probF pvalue32.30 ;
run;
proc export data = seq.de_sl_frag_&cell._FbyCase   
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex/de_sl_frag_&cell._FbyCase.csv"
dbms = csv replace ;
run;

%mend ;

%output (CD4) ;
%output (CD8) ;


