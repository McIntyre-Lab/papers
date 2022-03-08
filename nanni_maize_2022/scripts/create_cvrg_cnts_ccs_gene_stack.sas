

libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";




/* (1) create list of samples ccs data
        mo17 oz have 2 'rep's, rest only 1 'rep'

    make.list_ccs_samples


    (2) create stack coverage counts file for short read data

    make.cvrg_cnts_shrtRead
*/



filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2015/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/* list of samples to loop over - */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv"
out = pb_list2
dbms = csv replace ;
getnames = no ;
guessingrows = MAX ;
run;

data make.list_ccs_samples ;
set  pb_list2 ;
if var4 = "21-2" then rep = 2 ;
else rep = 1 ;
older = compress(var4||'_'||var2||'_'||var3) ;
oldSample = tranwrd(older, "-", "_") ;
sample = compress(var2||'_'||var3||'_'||rep) ;
rename var2 = geno ;
rename var3 = trt ;
keep sample oldSample var2 var3 rep;
run;

data pb_list ;
set make.list_ccs_samples ;
run ;

/* stack together */
%macro stacking (oldSample, sample, geno, trt, rep) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_ccs_b73_gene/cvrg_cnts_gene_&oldSample..csv"
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

%iterdataset(dataset=pb_list, function=%nrstr(%stacking(&oldSample, &sample, &geno, &trt, &rep);)); 


data make.cvrg_cnts_ccs ;
retain sampleID geno trt rep ;
length sampleID $ 12.;
length geno $ 5. ;
length rep $ 2. ;
set cnts2_: ;
run;


/* calculate mean of wt_tpm by genotype_trt  */
data ccs_4_means ;
retain sampleID geno_trt ;
set make.cvrg_cnts_ccs ;
geno_trt = compress(geno||'_'||trt) ;
run ;

proc sort data = ccs_4_means ;
by primary_FBgn geno_trt ;
run ;
 
proc means data = ccs_4_means  stackodsoutput  ;
by primary_FBgn geno_trt ;
var wt_tpm ;
id geno_trt ;
ods output summary = means_tpm_ccs ;
run ;

proc transpose data = means_tpm_ccs out = flip_means_tpm_ccs prefix = mean_tpm_ ;
by primary_FBgn ;
var mean ;
id geno_trt ;
run;

data means_tpm_ccs_sbys ;
set flip_means_tpm_ccs ;
drop _name_ ;
run ;

data make.means_cvrg_cnts_tpm_ccs_sbys;
set means_tpm_ccs_sbys;
run ;



/* sbys by genotype - 

    use 'old' name 
    drop bad samples  */
data pblist ;
set make.list_ccs_samples ;
keep oldSample sample ;
rename sample = sampleID ;
run ;

data pb_cc_shrt;
set make.cvrg_cnts_ccs ;
run;

proc sort data = pblist ;
by sampleID ;
proc sort data = pb_cc_shrt ;
by sampleID;
run;

data cc_ccs_rename ;
merge pblist (in=in2) pb_cc_shrt (in=in2) ;
by sampleID ;
run ;

data cc2_ccs_rename ;
retain oldSample sampleID ;
set cc_ccs_rename ;
drop sampleID ;
rename oldSample = sampleID ;
run ;


%macro looping (geno) ;

data cc_ccs_&geno. ;
set cc2_ccs_rename ;
where geno = "&geno." ;
run ;

proc sort data = cc_ccs_&geno. ;
by primary_FBgn ;
run ;

proc transpose data = cc_ccs_&geno. out = flip_ccs_&geno. suffix = _tpm ;
by primary_FBgn ;
var wt_tpm ;
id sampleID ;
run;

data flip2_ccs_&geno. ;
set flip_ccs_&geno. ;
drop _name_ ;
run ;

data cvrg_ccs_&geno._sbys ;
set flip2_ccs_&geno. ;
run ;

data means_&geno ;
set make.means_cvrg_cnts_tpm_ccs_sbys ;
keep primary_FBgn mean_tpm_&geno._Amb mean_tpm_&geno._Ele ;
run ;

data cvrg_ccs_tpm_&geno._sbys ;
merge means_&geno (in=in1) cvrg_ccs_&geno._sbys (in=in2) ;
by primary_FBgn ;
run ;

data make.cvrg_ccs_tpm_&geno._sbys ;
set cvrg_ccs_tpm_&geno._sbys ;
run;

data make.design_&geno._ccs ;
set cc_ccs_&geno. ;
rename trt = treatment ;
rename geno = genotype ;
keep sampleID trt geno ;
run ;

proc sort data = make.design_&geno._ccs nodups ;
by _all_ ;
run;


proc export data = make.cvrg_ccs_tpm_&geno._sbys 
outfile ="/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvrg_ccs_tpm_&geno._sbys.csv"
dbms = csv replace ;
run;

proc export data = make.design_&geno._ccs 
outfile ="/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/design_&geno._ccs.csv"
dbms = csv replace ;
run;

%mend ;

%looping (b73) ;
%looping (c123) ;
%looping (hp301) ;
%looping (mo17) ;
%looping (nc338) ;

