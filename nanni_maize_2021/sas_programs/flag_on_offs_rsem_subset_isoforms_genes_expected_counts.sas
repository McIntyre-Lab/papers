
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
 
/* on off calls --> EXPECTED COUNTS for Isoforms and Genes
            for cnt >0 and cnt > 5

    merge in design

chamber     trt
C1          amb
C2          amb
       
C4          ele
C5          ele

C6          amb
C7          ele

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


%macro on_offs (type, shrt, id, cnt) ;

proc sort data = pacbio.rsem_subset_&type._expCnt_stk;
by sample ;
proc sort data = design ;
by sample ;
run ;

data stacked_cnts_&type._&cnt oops ;
merge design (in=in1) pacbio.rsem_subset_&type._expCnt_stk (in=in2) ;
by sample ;
if in1 and in2 then output stacked_cnts_&type._&cnt ;
else output oops ;
run;

data stacked_&type._&cnt ;
retain &id.ID  sample genotype plant chamber trt expected_count ;
set stacked_cnts_&type._&cnt ;
run ;

proc sort data = stacked_&type._&cnt ;
by genotype trt plant chamber ;
run ;

data stacked_on_&type._&cnt ;
set stacked_&type._&cnt ;
if expected_count > 0 then &shrt._on0 = 1 ; else &shrt._on0 =0;
if expected_count > 5 then &shrt._on5 = 1 ; else &shrt._on5 =0;
run ;
/*
proc freq data = stacked_on_&type._&cnt;
by genotype ;
tables plant chamber * trt ;
run;
*/
proc sort data = stacked_on_&type._&cnt ;
by genotype trt &id.ID;
run ;

proc means data = stacked_on_&type._&cnt noprint ;
by genotype trt &id.ID ;
var &shrt._on&cnt. ;
output out = trt_on&cnt._&type. mean=trt_percent_on ;
run ;

proc sort data = trt_on&cnt._&type. ;
by &id.ID ;
run ;

proc transpose data = trt_on&cnt._&type. out = trt_on_sbys&cnt._&type. delimiter = _ ;
by &id.ID ;
id genotype trt ;
var trt_percent_on ;
run ;

data pacbio.sub_geno_trt_&shrt._onCall_Cnt&cnt. ;
set trt_on_sbys&cnt._&type. ;
by &id.ID;
/* B73 */
if B73_Amb > 0.5 then flag_B73_Amb_on&cnt. = 1 ; else flag_B73_Amb_on&cnt. = 0;
if B73_Ele > 0.5 then flag_B73_Ele_on&cnt. = 1 ; else flag_B73_Ele_on&cnt. = 0;

if flag_B73_Amb_on&cnt. = 0 and flag_B73_Ele_on&cnt. = 0 then flag_B73_trt_trnscpt_on&cnt. = 0; else flag_B73_trt_trnscpt_on&cnt. = 1;
if flag_B73_Amb_on&cnt. = 1 and flag_B73_Ele_on&cnt. = 1 then flag_B73_trt_trnscpt_all_on&cnt. = 1; else flag_B73_trt_trnscpt_all_on&cnt. = 0;

if B73_Amb = 1 then flag_B73_Amb_on&cnt._allReps = 1 ; else flag_B73_Amb_on&cnt._allReps = 0;
if B73_Ele = 1  then flag_B73_Ele_on&cnt._allReps = 1 ; else flag_B73_Ele_on&cnt._allReps = 0;

if B73_Amb = 0 then flag_B73_Amb_off&cnt._allReps = 1 ; else flag_B73_Amb_off&cnt._allReps = 0;
if B73_Ele = 0  then flag_B73_Ele_off&cnt._allReps = 1 ; else flag_B73_Ele_off&cnt._allReps = 0;
/* Mo17 */
if Mo17_Amb > 0.5 then flag_Mo17_Amb_on&cnt. = 1 ; else flag_Mo17_Amb_on&cnt. = 0;
if Mo17_Ele > 0.5 then flag_Mo17_Ele_on&cnt. = 1 ; else flag_Mo17_Ele_on&cnt. = 0;

if flag_Mo17_Amb_on&cnt. = 0 and flag_Mo17_Ele_on&cnt. = 0 then flag_Mo17_trt_trnscpt_on&cnt. = 0; else flag_Mo17_trt_trnscpt_on&cnt. = 1;
if flag_Mo17_Amb_on&cnt. = 1 and flag_Mo17_Ele_on&cnt. = 1 then flag_Mo17_trt_trnscpt_all_on&cnt. = 1; else flag_Mo17_trt_trnscpt_all_on&cnt. = 0;

if Mo17_Amb = 1 then flag_Mo17_Amb_on&cnt._allReps = 1 ; else flag_Mo17_Amb_on&cnt._allReps = 0;
if Mo17_Ele = 1  then flag_Mo17_Ele_on&cnt._allReps = 1 ; else flag_Mo17_Ele_on&cnt._allReps = 0;

if Mo17_Amb = 0 then flag_Mo17_Amb_off&cnt._allReps = 1 ; else flag_Mo17_Amb_off&cnt._allReps = 0;
if Mo17_Ele = 0  then flag_Mo17_Ele_off&cnt._allReps = 1 ; else flag_Mo17_Ele_off&cnt._allReps = 0;
/* C123 */
if C123_Amb > 0.5 then flag_C123_Amb_on&cnt. = 1 ; else flag_C123_Amb_on&cnt. = 0;
if C123_Ele > 0.5 then flag_C123_Ele_on&cnt. = 1 ; else flag_C123_Ele_on&cnt. = 0;

if flag_C123_Amb_on&cnt. = 0 and flag_C123_Ele_on&cnt. = 0 then flag_C123_trt_trnscpt_on&cnt. = 0; else flag_C123_trt_trnscpt_on&cnt. = 1;
if flag_C123_Amb_on&cnt. = 1 and flag_C123_Ele_on&cnt. = 1 then flag_C123_trt_trnscpt_all_on&cnt. = 1; else flag_C123_trt_trnscpt_all_on&cnt. = 0;

if C123_Amb = 1 then flag_C123_Amb_on&cnt._allReps = 1 ; else flag_C123_Amb_on&cnt._allReps = 0;
if C123_Ele = 1  then flag_C123_Ele_on&cnt._allReps = 1 ; else flag_C123_Ele_on&cnt._allReps = 0;

if C123_Amb = 0 then flag_C123_Amb_off&cnt._allReps = 1 ; else flag_C123_Amb_off&cnt._allReps = 0;
if C123_Ele = 0  then flag_C123_Ele_off&cnt._allReps = 1 ; else flag_C123_Ele_off&cnt._allReps = 0;
/* NC338 */
if NC338_Amb > 0.5 then flag_NC338_Amb_on&cnt. = 1 ; else flag_NC338_Amb_on&cnt. = 0;
if NC338_Ele > 0.5 then flag_NC338_Ele_on&cnt. = 1 ; else flag_NC338_Ele_on&cnt. = 0;

if flag_NC338_Amb_on&cnt. = 0 and flag_NC338_Ele_on&cnt. = 0 then flag_NC338_trt_trnscpt_on&cnt. = 0; else flag_NC338_trt_trnscpt_on&cnt. = 1;
if flag_NC338_Amb_on&cnt. = 1 and flag_NC338_Ele_on&cnt. = 1 then flag_NC338_trt_trnscpt_all_on&cnt. = 1; else flag_NC338_trt_trnscpt_all_on&cnt. = 0;

if NC338_Amb = 1 then flag_NC338_Amb_on&cnt._allReps = 1 ; else flag_NC338_Amb_on&cnt._allReps = 0;
if NC338_Ele = 1  then flag_NC338_Ele_on&cnt._allReps = 1 ; else flag_NC338_Ele_on&cnt._allReps = 0;

if NC338_Amb = 0 then flag_NC338_Amb_off&cnt._allReps = 1 ; else flag_NC338_Amb_off&cnt._allReps = 0;
if NC338_Ele = 0  then flag_NC338_Ele_off&cnt._allReps = 1 ; else flag_NC338_Ele_off&cnt._allReps = 0;
/* Hp301 */
if Hp301_Amb > 0.5 then flag_Hp301_Amb_on&cnt. = 1 ; else flag_Hp301_Amb_on&cnt. = 0;
if Hp301_Ele > 0.5 then flag_Hp301_Ele_on&cnt. = 1 ; else flag_Hp301_Ele_on&cnt. = 0;

if flag_Hp301_Amb_on&cnt. = 0 and flag_Hp301_Ele_on&cnt. = 0 then flag_Hp301_trt_trnscpt_on&cnt. = 0; else flag_Hp301_trt_trnscpt_on&cnt. = 1;
if flag_Hp301_Amb_on&cnt. = 1 and flag_Hp301_Ele_on&cnt. = 1 then flag_Hp301_trt_trnscpt_all_on&cnt. = 1; else flag_Hp301_trt_trnscpt_all_on&cnt. = 0;

if Hp301_Amb = 1 then flag_Hp301_Amb_on&cnt._allReps = 1 ; else flag_Hp301_Amb_on&cnt._allReps = 0;
if Hp301_Ele = 1  then flag_Hp301_Ele_on&cnt._allReps = 1 ; else flag_Hp301_Ele_on&cnt._allReps = 0;

if Hp301_Amb = 0 then flag_Hp301_Amb_off&cnt._allReps = 1 ; else flag_Hp301_Amb_off&cnt._allReps = 0;
if Hp301_Ele = 0  then flag_Hp301_Ele_off&cnt._allReps = 1 ; else flag_Hp301_Ele_off&cnt._allReps = 0;

keep &id.ID flag_: ;
run;
/*
proc freq data = pacbio.sub_geno_trt_&shrt._onCall_Cnt&cnt. ;
tables flag_: ;
run;
*/
%mend ;

%on_offs (isoforms, isoform, transcript, 0) ;
%on_offs (genes, gene, gene, 0) ;
%on_offs (isoforms, isoform, transcript, 5) ;
%on_offs (genes, gene, gene, 5) ;

/* export on off flags */

%macro exporting (shrt, tpm) ;

proc export data = pacbio.sub_geno_trt_&shrt._onCall_Cnt&cnt.
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/sub_geno_trt_&shrt._onCall_expCnt&cnt..tsv"
dbms = tab replace ;
run  ;
%mend ;

%exporting (isoform, 0);
%exporting (gene, 0);
%exporting (isoform, 5);
%exporting (gene, 5);






