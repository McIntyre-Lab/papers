

libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";




/* (1) create list of samples short read data
    convert plant # to rep # after sortingn by gene and trg
    flag bad samples

    make.list_shrtRd_samples_plant2rep


    (2) create stack coverage counts file for short read data
        where flag_bad_sample = 0 

    make.cvrg_cnts_shrtRead
*/


filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);


/* create list of samples to loop over  */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/maize_rnaseq_samples_fix_noHeader.csv"
out = list3
dbms = csv replace ;
getnames = no ;
guessingrows = MAX ;
run;

data list2 ;
set  list3 ;
geno = scan(var1, 1, '_') ;
trt = scan(var1, 4, '_') ;
plant = scan(var1, 2, '_') ; /* convert plant into arbitrary rep so easier to loop */
rename var1 = oldSample ;
run;

proc sort data = list2 nodups ;
by geno trt ;
run;

data make.list_shrtRd_samples_plant2rep ;
retain oldSample sample geno trt rep ;
set list2 ;
rep + 1 ;
by geno trt ;
if first.geno or first.trt then rep = 1 ;
sample = compress(geno||'_'||trt||'_'||rep) ; 
if oldSample = "B73_P4_C1_Amb" then flag_bad_sample = 1 ;
else if oldSample = "B73_P4_C6_Amb" then flag_bad_sample = 1 ;
else if oldSample = "B73_P1_C7_Ele" then flag_bad_sample = 1 ;
else flag_bad_sample = 0 ;

run;




/* list of samples to loop over - ex: mo17_p4_c6_amb */
data list ;
set make.list_shrtRd_samples_plant2rep ; 
where flag_bad_sample = 0 ;
drop flag_bad_sample ;
run ;


/* stack together */
%macro stacking (oldSample, sample, geno, trt, rep) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_shortRead_b73_gene/cvrg_cnts_gene_&oldSample..csv"
out = cnts_&sample.
dbms = csv replace ;
guessingrows = MAX ;
run;

data cnts2_&sample.;
retain sampleID ;
set cnts_&sample ;
sampleID = "&sample." ;
geno = "&geno." ;
trt = "&trt";
rep = "&rep." ;
run;

proc sort data = cnts2_&sample. ;
by primary_FBgn ;
run;

%mend ;
/*%stacking (B73_P1_C1_Amb, B73_Amb_1);*/

%iterdataset(dataset=list, function=%nrstr(%stacking(&oldSample, &sample, &geno, &trt, &rep);)); 


data make.cvrg_cnts_shrtRead ;
retain sampleID geno trt rep ;
length sampleID $ 12.;
length geno $ 5. ;
length rep $ 2. ;
set cnts2_: ;
run;


/* calculate mean of wt_tpm by genotype_trt  */
data shrt_4_means ;
retain sampleID geno_trt ;
set make.cvrg_cnts_shrtRead ;
geno_trt = compress(geno||'_'||trt) ;
run ;

proc sort data = shrt_4_means ;
by primary_FBgn geno_trt ;
run ;
 
proc means data = shrt_4_means  stackodsoutput  ;
by primary_FBgn geno_trt ;
var wt_tpm ;
id geno_trt ;
ods output summary = means_tpm_shrtRd ;
run ;

proc transpose data = means_tpm_shrtRd out = flip_means_tpm_shrtRd prefix = mean_tpm_ ;
by primary_FBgn ;
var mean ;
id geno_trt ;
run;

data means_tpm_shrtRd_sbys ;
set flip_means_tpm_shrtRd ;
drop _name_ ;
run ;

data make.means_cvrg_cnts_tpm_shrtRd_sbys;
set means_tpm_shrtRd_sbys;
run ;


/* sbys by genotype - 

    use 'old' name 
    drop bad samples  */
data list ;
set make.list_shrtRd_samples_plant2rep ;
keep oldSample sample flag_bad_sample;
rename sample = sampleID ;
run ;

data cc_shrt;
set make.cvrg_cnts_shrtRead ;
run;

proc sort data = list ;
by sampleID ;
proc sort data = cc_shrt ;
by sampleID;
run;

data cc_shrt_rename ;
merge list (in=in2) cc_shrt (in=in2) ;
by sampleID ;
run ;

data cc_renamed ;
retain oldSample sampleID ;
set cc_shrt_rename ;
where flag_bad_sample = 0 ;
drop sampleID ;
rename oldSample = sampleID ;
run ;


%macro looping (geno) ;

data cc_&geno. ;
set cc_renamed ;
where geno = "&geno." ;
run ;

proc sort data = cc_&geno. ;
by primary_FBgn ;
run ;

proc transpose data = cc_&geno. out = flip_&geno. suffix = _tpm ;
by primary_FBgn ;
var wt_tpm ;
id sampleID ;
run;

data flip2_&geno. ;
set flip_&geno. ;
drop _name_ ;
run ;

data cvrg_shrtRead_&geno._sbys ;
set flip2_&geno. ;
run ;

data means_&geno ;
set make.means_cvrg_cnts_tpm_shrtrd_sbys ;
keep primary_FBgn mean_tpm_&geno._Amb mean_tpm_&geno._Ele ;
run ;

data cvrg_shrtRead_tpm_&geno._sbys ;
merge means_&geno (in=in1) cvrg_shrtRead_&geno._sbys (in=in2) ;
by primary_FBgn ;
run ;

data make.cvrg_shrtRead_tpm_&geno._sbys ;
set cvrg_shrtRead_tpm_&geno._sbys ;
run;

data make.design_&geno._shrtRd ;
set cc_&geno. ;
plant = compress(scan(sampleID, 2, "_")) ;
chamber = compress(scan(sampleID, 3, "_")) ;
rename trt = treatment ;
rename geno = genotype ;
keep sampleID geno trt plant chamber ;
run ;

proc sort data = make.design_&geno._shrtRd nodups ;
by _all_ ;
run;

proc export data = make.cvrg_shrtRead_tpm_&geno._sbys
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvrg_shrtRead_tpm_&geno._sbys.csv"
dbms = csv replace ;
run;

proc export data = make.design_&geno._shrtRd
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/design_&geno._shrtRd.csv"
dbms = csv replace ;
run;

%mend ;

%looping (B73) ;
%looping (C123) ;
%looping (Hp301) ;
%looping (Mo17) ;
%looping (NC338) ;


