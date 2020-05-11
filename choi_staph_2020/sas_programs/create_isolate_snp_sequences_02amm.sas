libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
libname TB14 "/home/ammorse/TB14/staph_relapse/sasdata";

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);


%macro lists (ref) ;


data list_&ref ;
set relapse.sampleID_st_cc_ref ;
where reference = "&ref." ;
keep sampleID reference pairNum;
rename reference = ref ;
run ;
%mend ;

%lists (MSSA476) ;
%lists (ED98) ;
%lists (Newman) ;
%lists (ST20130941) ;
%lists (TCH60) ;
%lists (CA_347) ;

/* drop bad sample and pair from ED98 list*/

data list_ed98 ;
set list_ed98;
if sampleID = 3938 or sampleID = 3867 then delete ;
run ;
 
/* import vcf without header      
        create genotype calls for each isolate
        missing = -1
        ref = 0 */
%macro import_vcf (ref, sampleID, pairNum) ;

proc import datafile = "/home/ammorse/TB14/staph_relapse/vcf_w_targets_file/&sampleID._2_CCRef_&ref._renamed_snps_noHeader.vcf"
out = vcf_&sampleID._&ref. 
dbms = tab replace ;
guessingrows = MAX ;
run;

data &ref._&sampleID._vcf ;
set vcf_&sampleID._&ref;
if find(INFO, "INDEL") ge 1 then flag_indel = 1 ;
else flag_indel = 0 ;
rename _chrom = chrom ;
rename _&sampleID = isolate_&sampleID ;
keep _chrom pos ref alt _&sampleID flag_indel;
run ;

data seq_&ref._&sampleID. ;
length seq_&sampleID $4. ;
set &ref._&sampleID._vcf ;
where flag_indel = 0 ;
if isolate_&sampleID =: "1:" then seq_&sampleID  = alt ;
else if isolate_&sampleID =: "0:" then seq_&sampleID  = ref ;
else if isolate_&sampleID =: "." then seq_&sampleID  = "M" ;
    else seq_&sampleID = "M";
run;

proc sort data = seq_&ref._&sampleID. ;
by chrom pos ;
run ;

%mend ;

%iterdataset(dataset=list_mssa476, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));
%iterdataset(dataset=list_ed98, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));
%iterdataset(dataset=list_newman, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));
%iterdataset(dataset=list_ca_347, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));
%iterdataset(dataset=list_tch60, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));
%iterdataset(dataset=list_st20130941, function=%nrstr(%import_vcf(&ref, &sampleID, &pairNum);));

/* target file */
%macro targets (ref) ;

proc import datafile = "/home/ammorse/TB14/staph_relapse/vcf_by_ref/&ref._snps.tsv"
out = target2_&ref. 
dbms = tab replace ;
getnames = no;
run ;

data target_&ref ;
set  target2_&ref ;
rename var1 = chrom ;
rename var2 = pos ;
ref = scan(var3, 1, ',');
alt = scan(var3, 2, ',');
drop var3 ;
run ;

proc sort data = target_&ref. ;
by chrom pos ;
run;

%mend ;

%targets (MSSA476) ;
%targets (ED98) ;
%targets (Newman) ;
%targets (CA_347) ;
%targets (TCH60) ;
%targets (ST20130941) ;

/* calc number of obs (snps) for catting bases */
%macro numbers (ref) ;

data _NULL_ ;
if 0 then set target_&ref nobs = n ;
call symputx('nrows', n);
stop;
run;
%put ref=&ref;
%put nobs=&nrows;
%mend ;

%numbers (MSSA476);  /* 1105 */
%numbers (ED98) ;  /* 5746 */
%numbers (Newman) ;  /* 6725 */
%numbers (CA_347) ;  /* 1793 */
%numbers (TCH60) ;  /* 1950 */
%numbers (ST20130941) ;  /* 737 */


%macro add_target (ref, sampleID, pairNum, snpNum) ;

data add_&ref._&sampleID._seqs ;
merge seq_&ref._&sampleID. target_&ref. ;
by chrom pos ;
run ;

data add2_&ref._&sampleID._seqs ;
set  add_&ref._&sampleID._seqs ;
if seq_&sampleID = "" then seq_&sampleID = "M";
if flag_indel = 1 then delete ;
drop flag_indel;
run;

proc transpose data = add2_&ref._&sampleID._seqs out = flip_&ref._&sampleID.;
var seq_&sampleID ;
run;

data snp_&ref._&sampleID. ;
format fa $ &snpNum.. ;
informat  fa $ &snpNum..;
length fa $ &snpNum..;
set flip_&ref._&sampleID.;
fa = catt(of col1-col&snpNum) ;
length = &snpNum ;
rename _name_ = sampleID ;
keep fa _name_ length ;
run ;

%mend ;

%iterdataset(dataset=list_mssa476, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 1105);));
%iterdataset(dataset=list_ed98, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 5746);));
%iterdataset(dataset=list_newman, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 6725);));
%iterdataset(dataset=list_ca_347, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 1793);));
%iterdataset(dataset=list_tch60, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 1950);));
%iterdataset(dataset=list_st20130941, function=%nrstr(%add_target(&ref, &sampleID, &pairNum, 737);));

proc contents data = relapse.fasta_snps_mssa476 ; run;

%macro combine (ref, snpNum) ;

data combined_&ref. ;
set snp_&ref._: ;
run ;

data relapse.fasta_snps_&ref. ;
retain header fasta ;
format fasta $ &snpNum.. ;
informat fasta $ &snpNum..;
length fasta $ &snpNum..;
label header = "header" ;
set combined_&ref.;
ref = "&ref." ;
isolate = compress(scan(sampleID, 2, '_')) ;
header = compress(">snps_"||isolate||"_ref_"||ref||"_length_"||length) ;
fasta = tranwrd(fa, "M", "-");
keep header fasta ;
run ;


proc export data = relapse.fasta_snps_&ref. 
outfile = "!MCLAB/staph/seoung_ho_project/output/fasta_snps_&ref..csv"
dbms = csv replace ;
putnames = no ;
run;
%mend ;

%combine (MSSA476, 1105) ;
%combine (ED98, 5746) ;
%combine (Newman, 6725) ;
%combine (CA_347, 1793) ;
%combine (TCH60, 1950) ;
%combine (ST20130941, 737) ;



