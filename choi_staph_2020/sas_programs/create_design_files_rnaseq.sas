

libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";



/*
import list of fq files got from Batu
    80 samples but 1 failed sequencing so 79
    79 x 2 (PE) = 158 FQ files
    
Bacterial DNA Fowler-Seoung for WGS.xlsx 
    contains link between S No (sample number) - this is ID in FQ name - and sampleID
    sample number is the 1st '_' delimited value after the dash 
        e.g. the S58 in: 6124-S58_S83_L002_R2_001.fastq.gz (sample nums only go to 80 so cant be the s83)
        
Copy of Recurrent SAB isolates_for Lauren_200218.xlsx
    identifies the pairs of isolates (by sampleID)
           
*/


/* (1) list of fq files */
proc import datafile = "!MCLAB/staph/seoung_ho_project/design_files/fq_list.txt"
out = list
dbms = tab replace ;
getnames = no ;
guessingrows = MAX ;
run;

data relapse.fq_list_w_sampleNum ;
retain sampleNum var1 FQ read;
length sampleNum $5. read $3. ;
set list ;
rename var1 = fqName;
num1 = scan(var1, 2, '-');
sampleNum = compress(scan(num1, 1, '_'));
FQ = compress(scan(num1, 1, '.'));
read = compress(scan(num1, 4, '_'));
drop num1 ;
run ;  /* 158 = 79 samples x 2 (PE) */

/* (2) sampleNum to sampleID */
proc import datafile = "!MCLAB/staph/seoung_ho_project/Bacterial_DNA_Fowler-Seoung_for_WGS.tsv"
out = link
dbms = tab replace ;
guessingrows = MAX ;
run;

data sampleNum_sampleID ;
retain sampleNum sample_ID S_No ;
set link ;
sampleNum = compress("S"||S_NO);
rename sample_ID = sampleID ;
keep sampleNum S_No sample_ID ;
/* fix according to email on 24 Feb 2020 sample S44 should be 4553 */
if S_No = 44 then do;
    sampleNum = "S44";
    sample_ID = 4553 ;
end ;
run ;

proc sort data = relapse.fq_list_w_sampleNum ;
by sampleNum ;
proc sort data = sampleNum_sampleID  ;
by sampleNum ;
run ;

data find_bad ;
merge relapse.fq_list_w_sampleNum (in=in1) sampleNum_sampleID (in=in2) ;
by sampleNum ;
run ;
    /* S15 has no fqName, checked duke seq report and S15 NOT included --> this must be the failed library! */

data sampleNum_2_sampleID ;
set find_bad ;
if fqName = '' then flag_failed_library = 1 ;
else flag_failed_library = 0;
run ;

data sampleNum_2_sampleID_2 ;
set sampleNum_2_sampleID ;
drop fqName FQ read ;
run ;

proc sort data = sampleNum_2_sampleID_2 nodups ;
by _all_ ;
run;

data relapse.sampleNum_sampleID;
set sampleNum_2_sampleID_2 ;
run ;

proc freq data = relapse.sampleNum_sampleID ;
tables flag_failed_library ;
run;

/* fq list with sampleNum and sampleID */
data relapse.fq_list_w_sampleNum_sampleID ;
set sampleNum_2_sampleID ;
drop S_No ;
where flag_failed_library = 0;
run;  /* 158 obs!! */



/* (3) sampleID relapse pairs */
proc import datafile = "!MCLAB/staph/seoung_ho_project/Recurrent_SAB_Isolates_for_Lauren_200218_40Pairs.tsv"
out = pair
dbms = tab replace ;
guessingrows = MAX ;
run;

data relapse.isolate_pairs ;
set pair ;
rename interval__days_ = interval_in_days ; 
run;

data relapse.isolate_pairs_only ;
set  relapse.isolate_pairs ;
keep isolate_1 isolate_2 ;
run ;


/* export files */

%macro exporting (dataout) ;

proc export data = relapse.&dataout. 
outfile = "!MCLAB/staph/seoung_ho_project/design_files/&dataout..tsv"
dbms = tab replace ;
run ;
%mend ;

%exporting (sampleNum_sampleID);
%exporting (isolate_pairs);
%exporting (isolate_pairs_only);
%exporting (fq_list_w_sampleNum);
%exporting (fq_list_w_sampleNum_sampleID);

