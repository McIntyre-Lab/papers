libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;


/*


flagAnalyze_frag_DE = 1

model DE by cell type
  by fragment 
  t1d sex age  

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

/* combine CD4 and CD8 (excluding CD19) */

data ready_DE_control ;
set ready_DE_CD4 ready_DE_CD8 ;
where flag_case = 1;
run ;

proc sort data = ready_DE_control ;
by featureID ;
run;


proc glimmix data = ready_DE_control ;
by featureID ;
class flag_female population  ;
model log_uq_apn = enrollment_age flag_female population flag_female * population / htype =1 ;
lsmeans flag_female * population /slice = flag_female pdiff = all  ;
output out = seq.de_resid_frag_case resid=resid pred=pred student=student ;
ods output 
    tests1 = seq.de_t1_frag_case  
    lsmeans = seq.de_ls_frag_case 
    diffs = seq.de_pd_frag_case_FbyCell  
    slices = seq.de_sl_frag_case_FbyCell ;
    
run;

proc contents data = seq.de_sl_frag_case_FbyCell ; run;

%macro outputing (output, pvar) ;

proc datasets library = SEQ ;
modify &output. ;
format &pvar. pvalue32.30 ;
run;

proc export data = seq.&output.  
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_case_CD4_CD8/&output..csv"
dbms = csv replace ;
run;
%mend ;

%outputing (de_t1_frag_case, ProbF) ;
%outputing (de_ls_frag_case, Probt) ;
%outputing (de_pd_frag_case_FbyCell, Probt) ;
%outputing (de_sl_frag_case_FbyCell, ProbF ) ;


%macro outputing (output) ;

proc export data = seq.&output.  
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_case_CD4_CD8/&output..csv"
dbms = csv replace ;
run;
%mend ;
%outputing (de_resid_frag_case) ;
