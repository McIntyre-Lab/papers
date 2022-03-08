
/*
collapse to gene for (1) short read coverage counts


TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6

    Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
        rpk = region_depth / region_length
    Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
        scale_factor = (sum rpk for a sample) / 1,000,000
    Divide the RPK values by the “per million” scaling factor. This gives you TPM.
        tpm = rpk / scale_factor


collapse coverage counts to gene level
    files:  PROJ/rnaseq_cvrg_cnts_fusions

APN_gene = sum of all fusions * fusionLength divided by the sum of all fusionLengths
                       = (length1*apn1 + length2*apn2) / (length1+length2)

proc means use weight statement
	weight apn
	var length


*/

filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/* list fusions 2 genes */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_fusions.tsv"
out = anno 
dbms = tab replace ;
guessingrows = MAX ;
run; /* 402,691 rows  */


data fusion2gene ;
set anno ;
rename fusion_id = featureID ;
drop exon_id ;
where flag_multigene = 0 ;
run;  /* 389,033 features  */

proc sort data = fusion2gene dupout = chk nodups ;
by _all_ ;
run;  /* 193,394 fusions */


/* get fusion lengths from bed file */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_fusions_coverage.bed"
out = bed 
dbms = tab replace ;
getnames = no ;
guessingrows = MAX ;
run;  /* 196,394 fusions */

data fusionLen ;
set bed;
fusionLength = (var3 - var2 +1) ;
keep var4 fusionLength ;
rename var4 = featureID ;
run ;

/* add fusion lengths to fusion2gene */
proc sort data = fusion2gene ;
by featureID ;
proc sort data = fusionLen ;
by featureID ;
run;

data fusion2gene_len ;
merge fusion2gene (in=in1)  fusionLen (in=in2) ;
by featureID ;
if in1 and in2 ;
run ;  /* 193,394 fusions that are not multi-gene 
          count number genes as well */


/* list of samples to loop over - */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/maize_rnaseq_samples_fix_noHeader.csv"
out = list2
dbms = csv replace ;
getnames = no ;
guessingrows = MAX ;
run;

data list ;
set  list2 ;
rename var1 = sample ;
run;

proc sort data = list nodups ;
by sample;
run;

/* import coverage counts */
%macro summing (sample) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_shortRead_b73_fusions/cvrg_cnts_&sample..csv"
out = in_&sample.
dbms = csv replace ;
guessingrows = MAX;
run;   /* 196,059 fusions (has multigene)*/

proc sort data = fusion2gene_len ;
by featureID ;
proc sort data = in_&sample.;
by featureID ;
run;

data in2_&sample. ;
merge in_&sample. (in=in1) fusion2gene_len (in=in2) ;
by featureID ;
if in1 and in2;
run;  /* 193,394 fusions */


/* calculate TPM from rpkm 
    rpk = region_depth / region_length
    scale_factor = (sum rpk for a sample) / 1,000,000
    tpm = rpk / scale_factor  */

data in3_&sample. ;
set  in2_&sample. ;
rpk = region_depth / region_length ;
sample = "&sample.";
drop &sample. ;
run;

proc means data = in3_&sample. sum ;
var rpk;
output out = rpkSum_&sample. sum = ;
run;

proc sql noprint ;
select rpk into :rpk_sum
from rpkSum_&sample. ;
quit;

%put rpk is &rpk_sum ;

data in4_&sample. ;
retain sample featureID mapped_reads read_length region_length reads_in_region apn rpkm tpm mean std  ;
set in3_&sample. ;
scale_factor = (&rpk_sum / 1000000) ;
tpm = (rpk / scale_factor) ;
drop rpk scale_factor ;
run;


proc sort data = in4_&sample. ;
by primary_fbgn featureID ;
run;

/* output in stacked format - for read counts only */
proc means data = in4_&sample. sum stackods ;
by primary_fbgn;
ods output summary = want_&sample. ;
run ;

/* calc weighted ave for apn and tpm
	(length1*apn1 + length2*apn2) / (length1+length2)
*/
%macro weighting (var) ;

proc sql ;
create table wt_&var._&sample. as 
	select 
		primary_fbgn,
                sum(&var.*fusionLength) / sum(fusionLength) as wt_&var.
	from ( 
		select
		*
		from in4_&sample
		group by primary_fbgn)
	group by primary_fbgn;
	select * from wt_&var._&sample.;
quit ;


data wt2_&var._&sample.;
set wt_&var._&sample. ;
keep primary_fbgn wt_&var. ;
run ;

proc sort data = wt2_&var._&sample. ;
by primary_fbgn ;
run ;
%mend ;

%weighting (apn) ;
%weighting (tpm) ;


/* prep the count data */
data want2_&sample. ;
set want_&sample. ;
if find(variable, "apn") ge 1 then delete ;
if find(variable, "tpm") ge 1 then delete ;
if find(variable, "rpkm") ge 1 then delete ;
if find(variable, "mean") ge 1 then delete ;
if find(variable, "std") ge 1 then delete ;
if find(variable, "cv") ge 1 then delete ;
if find(variable, "flag_multigene") ge 1 then delete ;
if find(variable, "fusionLength") ge 1 then delete ;
run;

proc transpose data = want2_&sample. out = flip_&sample. ;
by primary_fbgn ;
id variable ;
run ;

data counts_&sample. ;
set flip_&sample. ;
drop _name_ ;
run ;

proc sort data = counts_&sample. ;
by primary_fbgn ;
run ;

/* merge count data and weighted APN values (already sorted) */
data counts_gene_&sample.  ;
merge counts_&sample. wt2_apn_&sample. wt2_tpm_&sample.;
by primary_fbgn ;
run ;

    data _null_;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    %let _EFIREC_ = 0;     /* clear export record count macro variable */
    file "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_shortRead_b73_gene/cvrg_cnts_gene_&sample..csv" 
    delimiter=',' DSD DROPOVER lrecl=32767;
    if _n_ = 1 then        /* write column names or labels */
     do;
       put
         "primary_FBgn" 
       ','
          "mapped_reads" 
       ','
          "read_length" 
       ','
          "region_length" 
       ','
          "reads_in_region" 
       ','
          "region_depth" 
       ','
          "wt_apn" 
       ','
          "wt_tpm" 
       ;
     end;
   set  COUNTS_GENE_&sample   end=EFIEOD;
       format primary_FBgn $15. ;
       format mapped_reads d12.3 ;
       format read_length d12.3 ;
       format region_length d12.3 ;
       format reads_in_region d12.3 ;
       format region_depth d12.3 ;
       format wt_apn best12. ;
       format wt_tpm best12. ;
     do;    
      EFIOUT + 1;
      put primary_FBgn $ @;
      put mapped_reads @;
      put read_length @;
      put region_length @;
      put reads_in_region @;
      put region_depth @;
      put wt_apn @;
      put wt_tpm ;
      ;
    end;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   if EFIEOD then call symputx('_EFIREC_',EFIOUT);
   run;

%mend ;

 /*%summing (Mo17_P4_C6_AMB) ; */


/*%iterdataset(dataset=test_list, function=%nrstr(%summing(&sample);)); */


%iterdataset(dataset=list, function=%nrstr(%summing(&sample);));






