libname sex "!MCLAB/sex_specific_splicing/sasdata";
libname local "/TB20/sex_specific_splicing/sasdata";

/* where is sas work 
proc sql;
   select setting 
   from sashelp.voption 
   where optname = 'WORK';
quit;
*/

%macro sumTR (species, genome) ;

/* sum tech reps for ujc counts */
data local.&species._cnts    ;
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "/TB20/sex_specific_splicing/&species._data_2_&genome._ujc_count.csv" delimiter = ','
MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat sampleID $26. ;
         informat geneID $26. ;
         informat jxnHash $256. ;
         informat numRead best32. ;
         format sampleID $26. ;
         format geneID $26. ;
         format jxnHash $256. ;
         format numRead best12. ;
      input
                  sampleID  $
                  geneID  $
                  jxnHash  $
                  numRead
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

data local.&species.Data2 ;
set local.&species._cnts;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
rep = scan(sampleID, 3, '_') ;
TR = scan(sampleID, 4, '_') ;
if find(sampleID, "dser") ge 1 then sample = compress(species||'_'||sex||'_'||rep) ;
else sample = compress("d"||species||'_'||sex||'_'||rep) ;
run ;

data sex.&species._jxnhash2geneID;
set  local.&species.Data2 ;
keep jxnHash geneID ;
run;

proc sort data = sex.&species._jxnhash2geneID nodups ;
by _all_ ;
run;

/* sum the tech reps  */
proc sort data = local.&species.Data2 ;
by sample jxnHash ;
run;

proc means data = local.&species.Data2 ;
by sample  jxnHash ;
var numRead ;
output out = local.&species._summed sum = sum ;
run;

data local.&species._summed2 ;
set local.&species._summed ;
drop _type_ _freq_ ;
rename sum = jxnHash_sumReadNum ;
run ;

proc sort data = local.&species._summed2 ;
by jxnHash ;
run;

/* make perm */
data sex.&species._ujc_cnts_sumTR_stack;
set local.&species._summed2 ;
run ;
%mend ;

%sumTR (dmel, dmel6) ;
%sumTR (dsim, dsim2) ;
%sumTR (dser, dser1) ;


%macro sumTR (species, genome) ;

/* sum tech reps for ujc counts */

data local.&species._cnts    ;
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "/TB20/sex_specific_splicing/&species._data_2_&genome._ujc_count.csv" delimiter = ','
MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat sampleID $26. ;
         informat geneID $26. ;
         informat jxnHash $256. ;
         informat numRead best32. ;
         format sampleID $26. ;
         format geneID $26. ;
         format jxnHash $256. ;
         format numRead best12. ;
      input
                  sampleID  $
                  geneID  $
                  jxnHash  $
                  numRead
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

data local.&species.Data2 ;
set local.&species._cnts;
species = scan(sampleID, 1, '_') ;
sex = scan(sampleID, 2, '_') ;
TR = scan(sampleID, 3, '_') ;
sample = compress(species||'_'||sex||'_rep1') ;
run ;

data sex.&species._jxnhash2geneID;
set  local.&species.Data2 ;
keep jxnHash geneID ;
run;

proc sort data = sex.&species._jxnhash2geneID nodups ;
by _all_ ;
run;

/* sum the tech reps  */
proc sort data = local.&species.Data2 ;
by sample geneID jxnHash ;
run;

proc means data = local.&species.Data2 ;
by sample geneID jxnHash ;
var numRead ;
output out = local.&species._summed sum = sum ;
run;

data local.&species._summed2 ;
set local.&species._summed ;
drop _type_ _freq_ ;
rename sum = jxnHash_sumReadNum ;
run ;

proc sort data = local.&species._summed2 ;
by jxnHash ;
run;

/* make perm */
data sex.&species._ujc_cnts_sumTR_stack;
set local.&species._summed2 ;
run ;
%mend ;

%sumTR (dsan, dsan1) ;
%sumTR (dyak, dyak2) ;


/* check above stack files are uniq on jxnHash */
ods pdf file="!MCLAB/sex_specific_splicing/check_sas_import_ujc_cnt_file_uniq_jxnHash.pdf";

%macro uniq_check (species) ;

proc sort data = local.&species.Data2 ;
by sampleID ;
run;

proc freq data = local.&species.Data2 noprint; 
by sampleID ;
tables jxnHash / out = check_&species. ;
run ;

data check2_&species. ;
set check_&species. ;
where count ne 1 ;
run;

/* check if all are 0 */
proc sql noprint;
    select count(*) into :nobs from check2_&species.;
quit;

/* If dataset is not empty, print the data */
%if &nobs > 1 %then %do;
    proc print data=check2_&species.;
        title "non-uniq jxnHash (count greater to 1) for &species. by sampleID ";
    run;
%end;
%else %do;
    data _null_;
        file print;
        put "jxnHash is uniq for &species. by sampleID in *_ujc_counts.csv input file";
    run;
%end;

%mend ;

%uniq_check (dmel) ;
%uniq_check (dsim) ;
%uniq_check (dser) ;
%uniq_check (dyak) ;
%uniq_check (dsan) ;

ods pdf close ;




%macro exping (species) ;

proc export data = sex.&species._ujc_cnts_sumTR_stack 
oufile = "!MCLAB/sex_specific_splicing/&species._ujc_cnts_sumTR_stack.csv"
dbms = csv replace ;
run ;
%mend ;

%exping (dmel) ;
%exping (dsim) ;
%exping (dser) ;
%exping (dsan) ;
%exping (dyak) ;



