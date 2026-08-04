libname sex "!MCLAB/sex_specific_splicing/sasdata";
libname local "/TB20/sex_specific_splicing/sasdata";

/* where is sas work 
proc sql;
   select setting 
   from sashelp.voption 
   where optname = 'WORK';
quit;
*/


%macro sumTR2 (species, genome) ;
  
   data LOCAL.&species._ERP    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile  "!MCLAB/sex_specific_splicing/erp_and_esp_output/fiveSpecies_2_&genome._ujc_er_vs_&species._data_2_&genome._ujc_noMultiGene_read_per_ERP.csv"
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sampleID $26. ;
        informat geneID $26. ;
        informat ERP $120. ;
        informat ERP_plus $125. ;
        informat flagDataOnlyExon best32. ;
        informat flagIR best32. ;
        informat numRead best32. ;
        format sampleID $26. ;
        format geneID $26. ;
        format ERP $120. ;
        format ERP_plus $125. ;
        format flagDataOnlyExon best12. ;
        format flagIR best12. ;
        format numRead best12. ;
     input
                 sampleID  $
                 geneID  $
                 ERP  $
                 ERP_plus  $
                 flagDataOnlyExon
                 flagIR
                 numRead
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data local.&species._erp2 ;
set local.&species._erp;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
rep = scan(sampleID, 3, '_') ;
TR = scan(sampleID, 4, '_') ;
if find(sampleID, "dser") ge 1 then sample = compress(species||'_'||sex||'_'||rep) ;
else sample = compress("d"||species||'_'||sex||'_'||rep) ;
run ;

/* sum the tech reps  */
proc sort data = local.&species._erp2 ;
by sample geneID ERP_plus;
run;

proc means data = local.&species._erp2 ;
by sample geneID ERP_plus;
var numRead ;
output out = local.&species._erp_sum sum = sum ;
run;

data local.&species._erp_sum2 ;
set local.&species._erp_sum ;
drop _type_ _freq_ ;
rename sum = ERP_sumReadNum ;
run ;

/* make perm */
data sex.&species._erp_cnts_sumTR_stack;
set local.&species._erp_sum2 ;
run ;


%mend ;

%sumTR2 (dmel, dmel6) ;
%sumTR2 (dsim, dsim2) ;
%sumTR2 (dser, dser1) ;


%macro sumTR2 (species, genome) ;

   data LOCAL.&species._ERP    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile  "!MCLAB/sex_specific_splicing/erp_and_esp_output/fiveSpecies_2_&genome._ujc_er_vs_&species._data_2_&genome._ujc_noMultiGene_read_per_ERP.csv"
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sampleID $26. ;
        informat geneID $26. ;
        informat ERP $120. ;
        informat ERP_plus $125. ;
        informat flagDataOnlyExon best32. ;
        informat flagIR best32. ;
        informat numRead best32. ;
        format sampleID $26. ;
        format geneID $26. ;
        format ERP $120. ;
        format ERP_plus $125. ;
        format flagDataOnlyExon best12. ;
        format flagIR best12. ;
        format numRead best12. ;
     input
                 sampleID  $
                 geneID  $
                 ERP  $
                 ERP_plus  $
                 flagDataOnlyExon
                 flagIR
                 numRead
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data local.&species._erp2 ;
set local.&species._erp;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
TR = scan(sampleID, 3, '_') ;
sample = compress(species||'_'||sex||'_rep1') ;
run ;

/* sum the tech reps  */
proc sort data = local.&species._erp2 ;
by sample geneID ERP_plus;
run;

proc means data = local.&species._erp2 ;
by sample geneID ERP_plus;
var numRead ;
output out = local.&species._erp_sum sum = sum ;
run;

data local.&species._erp_sum2 ;
set local.&species._erp_sum ;
drop _type_ _freq_ ;
rename sum = ERP_sumReadNum ;
run ;

/* make perm */
data sex.&species._erp_cnts_sumTR_stack;
set local.&species._erp_sum2 ;
run ;
%mend ;

%sumTR2 (dsan, dsan1) ;
%sumTR2 (dyak, dyak2) ;


%macro exping (species) ;

proc export data = sex.&species._erp_cnts_sumTR_stack 
oufile = "!MCLAB/sex_specific_splicing/&species._erp_cnts_sumTR_stack.csv"
dbms = csv replace ;
run ;
%mend ;

%exping (dmel) ;
%exping (dsim) ;
%exping (dser) ;
%exping (dsan) ;
%exping (dyak) ;






