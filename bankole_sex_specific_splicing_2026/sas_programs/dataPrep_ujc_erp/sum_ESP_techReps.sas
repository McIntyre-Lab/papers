libname sex "!MCLAB/sex_specific_splicing/sasdata";
libname local "/TB20/sex_specific_splicing/sasdata";

/* where is sas work 
proc sql;
   select setting 
   from sashelp.voption 
   where optname = 'WORK';
quit;
*/


%macro sumTR3 (species, genome) ;

     data LOCAL.&species._ESP    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/TB20/sex_specific_splicing/fiveSpecies_2_&genome._ujc_es_vs_&species._data_2_&genome._ujc_noMultiGene_read_per_ESP.csv"
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sampleID $26. ;
        informat ESP $163. ;
        informat flagDataOnlyExon best32. ;
        informat geneID $26. ;
        informat strand $1. ;
        informat seqname $26. ;
        informat numRead best32. ;
        format sampleID $26. ;
        format ESP $163. ;
        format flagDataOnlyExon best12. ;
        format geneID $26. ;
        format strand $1. ;
        format seqname $26. ;
        format numRead best12. ;
     input
                 sampleID  $
                 ESP  $
                 flagDataOnlyExon
                 geneID  $
                 strand  $
                 seqname  $
                 numRead
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;
/* sum tech reps for ujc counts */
proc import datafile = "/TB20/sex_specific_splicing/fiveSpecies_2_&genome._ujc_es_vs_&species._data_2_&genome._ujc_noMultiGene_read_per_ESP.csv"
out = local.&species._esp 
dbms = csv replace ;
guessingrows = MAX ;
run;

data local.&species._esp2 ;
set local.&species._esp;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
rep = scan(sampleID, 3, '_') ;
TR = scan(sampleID, 4, '_') ;
if find(sampleID, "dser") ge 1 then sample = compress(species||'_'||sex||'_'||rep) ;
else sample = compress("d"||species||'_'||sex||'_'||rep) ;
run ;

/* sum the tech reps  */
proc sort data = local.&species._esp2 ;
by sample geneID ESP flagDataOnlyExon;
run;

proc means data = local.&species._esp2 ;
by sample geneID ESP  flagDataOnlyExon;
var numRead ;
output out = local.&species._esp_sum sum = sum ;
run;

data local.&species._esp_sum2 ;
set local.&species._esp_sum ;
drop _type_ _freq_ ;
rename sum = ESP_sumReadNum ;
run ;

/* make perm */
data sex.&species._esp_cnts_sumTR_stack;
retain new ;
set local.&species._esp_sum2 ;
run ;


%mend ;

%sumTR3 (dmel, dmel6) ;
%sumTR3 (dsim, dsim2) ;
%sumTR3 (dser, dser1) ;


%macro sumTR3 (species, genome) ;
    data local.&species._esp    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/TB20/sex_specific_splicing/fiveSpecies_2_&genome._ujc_es_vs_&species._data_2_&genome._ujc_noMultiGene_read_per_ESP.csv"
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sampleID $26. ;
        informat ESP $163. ;
        informat flagDataOnlyExon best32. ;
        informat geneID $26. ;
        informat strand $1. ;
        informat seqname $26. ;
        informat numRead best32. ;
        format sampleID $26. ;
        format ESP $163. ;
        format flagDataOnlyExon best12. ;
        format geneID $26. ;
        format strand $1. ;
        format seqname $26. ;
        format numRead best12. ;
     input
                 sampleID  $
                 ESP  $
                 flagDataOnlyExon
                 geneID  $
                 strand  $
                 seqname  $
                 numRead
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data local.&species._esp2 ;
set local.&species._esp;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
TR = scan(sampleID, 3, '_') ;
sample = compress(species||'_'||sex||'_rep1') ;
run ;

/* sum the tech reps  */
proc sort data = local.&species._esp2 ;
by sample geneID ESP flagDataOnlyExon;
run;

proc means data = local.&species._esp2 ;
by sample geneID ESP flagDataOnlyExon;
var numRead ;
output out = local.&species._esp_sum sum = sum ;
run;

data local.&species._esp_sum2 ;
set local.&species._esp_sum ;
drop _type_ _freq_ ;
rename sum = ESP_sumReadNum ;
run ;

/* make perm */
data sex.&species._esp_cnts_sumTR_stack;
retain new ;
set local.&species._esp_sum2 ;
run ;

%mend ;


%sumTR3 (dsan, dsan1) ;
%sumTR3 (dyak, dyak2) ;


%macro exping (species) ;

proc export data = sex.&species._esp_cnts_sumTR_stack 
oufile = "!MCLAB/sex_specific_splicing/&species._esp_cnts_sumTR_stack.csv"
dbms = csv replace ;
run ;
%mend ;

%exping (dmel) ;
%exping (dsim) ;
%exping (dser) ;
%exping (dsan) ;
%exping (dyak) ;

