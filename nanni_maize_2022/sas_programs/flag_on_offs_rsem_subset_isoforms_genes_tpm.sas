
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
 
/* on off calls --> TPM for Isoforms and Genes
            for TPM >0 and TPM > 5

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


%macro on_offs (type, shrt, id, tpm) ;

proc sort data = pacbio.rsem_subset_&type._TPM_stk;
by sample ;
proc sort data = design ;
by sample ;
run ;

data stacked_cnts_&type._&tpm oops ;
merge design (in=in1) pacbio.rsem_subset_&type._TPM_stk (in=in2) ;
by sample ;
if in1 and in2 then output stacked_cnts_&type._&tpm ;
else output oops ;
run;

data stacked_&type._&tpm ;
retain &id.ID  sample genotype plant chamber trt tpm ;
set stacked_cnts_&type._&tpm ;
run ;

proc sort data = stacked_&type._&tpm ;
by genotype trt plant chamber ;
run ;

data stacked_on_&type._&tpm ;
set stacked_&type._&tpm ;
if tpm > 0 then &shrt._on0 = 1 ; else &shrt._on0 =0;
if tpm > 5 then &shrt._on5 = 1 ; else &shrt._on5 =0;
run ;
/*
proc freq data = stacked_on_&type._&tpm ;
by genotype ;
tables plant chamber * trt ;
run;
*/
proc sort data = stacked_on_&type._&tpm ;
by genotype trt &id.ID;
run ;

proc means data = stacked_on_&type._&tpm noprint ;
by genotype trt &id.ID ;
var &shrt._on&tpm. ;
output out = trt_on&tpm._&type. mean=trt_percent_on ;
run ;

proc sort data = trt_on&tpm._&type. ;
by &id.ID ;
run ;

proc transpose data = trt_on&tpm._&type. out = trt_on_sbys&tpm._&type. delimiter = _ ;
by &id.ID ;
id genotype trt ;
var trt_percent_on ;
run ;

data pacbio.sub_geno_trt_&shrt._onCall_tpm&tpm. ;
set trt_on_sbys&tpm._&type. ;
by &id.ID;
/* B73 */
if B73_Amb > 0.5 then flag_B73_Amb_on&tpm. = 1 ; else flag_B73_Amb_on&tpm. = 0;
if B73_Ele > 0.5 then flag_B73_Ele_on&tpm. = 1 ; else flag_B73_Ele_on&tpm. = 0;

if flag_B73_Amb_on&tpm. = 0 and flag_B73_Ele_on&tpm. = 0 then flag_B73_trt_trnscpt_on&tpm. = 0; else flag_B73_trt_trnscpt_on&tpm. = 1;
if flag_B73_Amb_on&tpm. = 1 and flag_B73_Ele_on&tpm. = 1 then flag_B73_trt_trnscpt_all_on&tpm. = 1; else flag_B73_trt_trnscpt_all_on&tpm. = 0;

if B73_Amb = 1 then flag_B73_Amb_on&tpm._allReps = 1 ; else flag_B73_Amb_on&tpm._allReps = 0;
if B73_Ele = 1 then flag_B73_Ele_on&tpm._allReps = 1 ; else flag_B73_Ele_on&tpm._allReps = 0;

if B73_Amb = 0 then flag_B73_Amb_off&tpm._allReps = 1 ; else flag_B73_Amb_off&tpm._allReps = 0;
if B73_Ele = 0 then flag_B73_Ele_off&tpm._allReps = 1 ; else flag_B73_Ele_off&tpm._allReps = 0;
/* Mo17 */
if Mo17_Amb > 0.5 then flag_Mo17_Amb_on&tpm. = 1 ; else flag_Mo17_Amb_on&tpm. = 0;
if Mo17_Ele > 0.5 then flag_Mo17_Ele_on&tpm. = 1 ; else flag_Mo17_Ele_on&tpm. = 0;

if flag_Mo17_Amb_on&tpm. = 0 and flag_Mo17_Ele_on&tpm. = 0 then flag_Mo17_trt_trnscpt_on&tpm. = 0; else flag_Mo17_trt_trnscpt_on&tpm. = 1;
if flag_Mo17_Amb_on&tpm. = 1 and flag_Mo17_Ele_on&tpm. = 1 then flag_Mo17_trt_trnscpt_all_on&tpm. = 1; else flag_Mo17_trt_trnscpt_all_on&tpm. = 0;

if Mo17_Amb = 1 then flag_Mo17_Amb_on&tpm._allReps = 1 ; else flag_Mo17_Amb_on&tpm._allReps = 0;
if Mo17_Ele = 1 then flag_Mo17_Ele_on&tpm._allReps = 1 ; else flag_Mo17_Ele_on&tpm._allReps = 0;

if Mo17_Amb = 0 then flag_Mo17_Amb_off&tpm._allReps = 1 ; else flag_Mo17_Amb_off&tpm._allReps = 0;
if Mo17_Ele = 0 then flag_Mo17_Ele_off&tpm._allReps = 1 ; else flag_Mo17_Ele_off&tpm._allReps = 0;
/* C123 */
if C123_Amb > 0.5 then flag_C123_Amb_on&tpm. = 1 ; else flag_C123_Amb_on&tpm. = 0;
if C123_Ele > 0.5 then flag_C123_Ele_on&tpm. = 1 ; else flag_C123_Ele_on&tpm. = 0;

if flag_C123_Amb_on&tpm. = 0 and flag_C123_Ele_on&tpm. = 0 then flag_C123_trt_trnscpt_on&tpm. = 0; else flag_C123_trt_trnscpt_on&tpm. = 1;
if flag_C123_Amb_on&tpm. = 1 and flag_C123_Ele_on&tpm. = 1 then flag_C123_trt_trnscpt_all_on&tpm. = 1; else flag_C123_trt_trnscpt_all_on&tpm. = 0;

if C123_Amb = 1 then flag_C123_Amb_on&tpm._allReps = 1 ; else flag_C123_Amb_on&tpm._allReps = 0;
if C123_Ele = 1  then flag_C123_Ele_on&tpm._allReps = 1 ; else flag_C123_Ele_on&tpm._allReps = 0;

if C123_Amb = 0 then flag_C123_Amb_off&tpm._allReps = 1 ; else flag_C123_Amb_off&tpm._allReps = 0;
if C123_Ele = 0 then flag_C123_Ele_off&tpm._allReps = 1 ; else flag_C123_Ele_off&tpm._allReps = 0;
/* NC338 */
if NC338_Amb > 0.5 then flag_NC338_Amb_on&tpm. = 1 ; else flag_NC338_Amb_on&tpm. = 0;
if NC338_Ele > 0.5 then flag_NC338_Ele_on&tpm. = 1 ; else flag_NC338_Ele_on&tpm. = 0;

if flag_NC338_Amb_on&tpm. = 0 and flag_NC338_Ele_on&tpm. = 0 then flag_NC338_trt_trnscpt_on&tpm. = 0; else flag_NC338_trt_trnscpt_on&tpm. = 1;
if flag_NC338_Amb_on&tpm. = 1 and flag_NC338_Ele_on&tpm. = 1 then flag_NC338_trt_trnscpt_all_on&tpm. = 1; else flag_NC338_trt_trnscpt_all_on&tpm. = 0;

if NC338_Amb = 1 then flag_NC338_Amb_on&tpm._allReps = 1 ; else flag_NC338_Amb_on&tpm._allReps = 0;
if NC338_Ele = 1 then flag_NC338_Ele_on&tpm._allReps = 1 ; else flag_NC338_Ele_on&tpm._allReps = 0;

if NC338_Amb = 0 then flag_NC338_Amb_off&tpm._allReps = 1 ; else flag_NC338_Amb_off&tpm._allReps = 0;
if NC338_Ele = 0 then flag_NC338_Ele_off&tpm._allReps = 1 ; else flag_NC338_Ele_off&tpm._allReps = 0;
/* Hp301 */
if Hp301_Amb > 0.5 then flag_Hp301_Amb_on&tpm. = 1 ; else flag_Hp301_Amb_on&tpm. = 0;
if Hp301_Ele > 0.5 then flag_Hp301_Ele_on&tpm. = 1 ; else flag_Hp301_Ele_on&tpm. = 0;

if flag_Hp301_Amb_on&tpm. = 0 and flag_Hp301_Ele_on&tpm. = 0 then flag_Hp301_trt_trnscpt_on&tpm. = 0; else flag_Hp301_trt_trnscpt_on&tpm. = 1;
if flag_Hp301_Amb_on&tpm. = 1 and flag_Hp301_Ele_on&tpm. = 1 then flag_Hp301_trt_trnscpt_all_on&tpm. = 1; else flag_Hp301_trt_trnscpt_all_on&tpm. = 0;

if Hp301_Amb = 1 then flag_Hp301_Amb_on&tpm._allReps = 1 ; else flag_Hp301_Amb_on&tpm._allReps = 0;
if Hp301_Ele = 1  then flag_Hp301_Ele_on&tpm._allReps = 1 ; else flag_Hp301_Ele_on&tpm._allReps = 0;

if Hp301_Amb = 0 then flag_Hp301_Amb_off&tpm._allReps = 1 ; else flag_Hp301_Amb_off&tpm._allReps = 0;
if Hp301_Ele = 0  then flag_Hp301_Ele_off&tpm._allReps = 1 ; else flag_Hp301_Ele_off&tpm._allReps = 0;

keep &id.ID flag_: ;
run;
/*
proc freq data = pacbio.sub_geno_trt_&shrt._onCall_tpm&tpm. ;
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

proc export data = pacbio.sub_geno_trt_&shrt._onCall_tpm&tpm.
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/sub_geno_trt_&shrt._onCall_tpm&tpm..tsv"
dbms = tab replace ;
run  ;
%mend ;

%exporting (isoform, 0);
%exporting (gene, 0);
%exporting (isoform, 5);
%exporting (gene, 5);



proc freq data =  pacbio.sub_geno_trt_gene_onCall_tpm0;
tables geneID / out = cnts_calls ;
run ;

data check_calls ;
set cnts_calls ;
where count ne  1 ;
run;



