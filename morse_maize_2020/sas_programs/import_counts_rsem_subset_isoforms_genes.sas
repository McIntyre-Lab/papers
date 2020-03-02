
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";



/* import counts --> FSM, ISM, NIC and NNC isoforms */

%macro importing (type, id, cnts, Cshrt) ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/combined_expression_&type._matrix.filtered_GTF.&cnts..txt"
out = &type._&cnts
dbms = tab replace ;
guessingrows = MAX ;
run;

/* 
create stack 
*/

proc transpose data = &type._&cnts out = &type._&cnts._stack ;
by &id._id ;
run;

data &type._&cnts._stack2 ;
set &type._&cnts._stack ;
rename col1 = &cnts;
rename _name_ = sample ;
rename &id._id = &id.ID ;
run ;

data pacbio.rsem_subset_&type._&Cshrt._stk ;
set  &type._&cnts._stack2 ;
label sample = "sample" ;
run ;

%mend ;

%importing (isoforms, transcript, TPM, tpm) ;
%importing (isoforms, transcript, expected_count, expCnt) ;

%importing (genes, gene, TPM, tpm) ;
%importing (genes, gene, expected_count, expCnt) ;

proc freq data = genes_tpm ;
tables gene_id / out = cnts ;
run;
data check ;
set cnts ;
where count ne 1 ;
run;  /* uniq */



