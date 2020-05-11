


libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
libname TB14 "/home/ammorse/TB14/staph_relapse/sasdata";

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);




proc import datafile = "/home/ammorse/TB14/staph_relapse/design_files/isolate_pairs_only_noHeader.csv"
out = list 
dbms = csv replace ;
getnames = no ;
run;

data list_pairs ;
set list ;
rename var1 = pair1 ;
rename var2 = pair2 ;
if var1 = 3938 or var2 = 3938 then delete ;  /* drop pair with failed library */
run;


%macro tables (pair1, pair2, depth) ;

proc import datafile = "/home/ammorse/TB14/staph_relapse/vcf_pairs/&pair1._to_&pair2._allele_frq_minDepth&depth..frq"
out = frq_&depth._&pair1._to_&pair2
dbms = tab replace ;
guessingrows = MAX ;
run;

data frq_DP&depth._&pair1._to_&pair2;
set frq_&depth._&pair1._to_&pair2 ;
allele_freq_refBase = 1*(scan(_allele_freq_, 2, ':'));
allele_freq_altBase1 = 1*(scan(var6, 2, ':'));
allele_freq_altBase2 = 1*(scan(var7, 2, ':'));
if allele_freq_altBase2 = '' then allele_freq_altBase2 = 0;
drop _allele_freq_ var6 var7 ;
run;

data frq2_DP&depth._&pair1._to_&pair2;
set frq_DP&depth._&pair1._to_&pair2 ;
flag_altBase_mult = 0;

if allele_freq_refBase = 1 then flag_refBase = 1 ;
    else flag_refBase = 0;
if allele_freq_altBase1 =1 then flag_altBase = 1 ;
else if (allele_freq_altBase1 + allele_freq_altBase2) =1 then do ;
    flag_altBase = 1 ;
    flag_altBase_mult = 1;
    end ;
else flag_altBase = 0 ; 
run;

proc freq data = frq2_DP&depth._&pair1._to_&pair2 ;
tables flag_refBase   / out = cnting_ref;
run;

proc freq data = frq2_DP&depth._&pair1._to_&pair2 ;
tables flag_altBase   / out = cnting_alt;
run;

data cnting_ref2 ;
length isolate_pair $28.;
label count = "count_flag_refBase";
label percent = "percent_flag_refBase" ;
set cnting_ref ;
rename count = count_flag_refBase ;
rename percent = percent_flag_refBase ;
mean_minimum_depth = &depth.   ;
isolate_pair = compress(('isolate_'||&pair1.||'_'||'isolate_'||&pair2.) );
run ;

data cnting_alt2 ;
length isolate_pair $28.;
label count = "count_flag_altBase";
label percent = "percent_flag_altBase" ;
set cnting_alt ;
rename count = count_flag_altBase ;
rename percent = percent_flag_altBase ;
mean_minimum_depth = &depth.   ;
isolate_pair = compress(('isolate_'||&pair1.||'_'||'isolate_'||&pair2.) );
run ;

proc freq data = frq2_DP&depth._&pair1._to_&pair2 ;
tables  flag_altBase_mult / out = cnting_mult;
run;

data cnting_mult2 ;
length isolate_pair $28.;
label count = "count_flag_altBase_mult" ;
set cnting_mult ;
isolate_pair = compress(('isolate_'||&pair1.||'_'||'isolate_'||&pair2.) );
if flag_altBase_mult = 1 ;
rename count = count_flag_altBase_mult ;
drop percent ;
run ;

data all_frq_minDP&depth._&pair1._to_&pair2 ;
merge cnting_ref2 cnting_alt2 cnting_mult2  ;
run;


data allele_frq_minDP&depth._&pair1._to_&pair2 ;
retain isolate_pair mean_minimum_depth total_snp_num percent_altBase mult_altBase_num ;
set all_frq_minDP&depth._&pair1._to_&pair2  ;
if count_flag_refBase = count_flag_altBase then total_snp_num = count_flag_refBase ;
else total_snp_num = -1 ;

if percent_flag_altBase = 100 and flag_altBase = 1 then percent_altBase = 100  ;
    else percent_altBase = percent_flag_altBase ;

if flag_altBase_mult ge 1 then mult_altBase_num = count_flag_altBase_mult ;
    else count_flag_altBase_mult =0 ;

drop count_flag_refBase 
    count_flag_altBase 
    percent_flag_refBase 
    percent_flag_altBase 
    flag_altBase 
    flag_altBase_mult 
    count_flag_altBase_mult flag_refBase ;
run;


%mend ;

/*  %tables (3978, 4020, 5) ;  */

%iterdataset(dataset=list_pairs, function=%nrstr(%tables(&pair1, &pair2, 5);));
%iterdataset(dataset=list_pairs, function=%nrstr(%tables(&pair1, &pair2, 10);));


proc sql  ; 
select memname into :deletelist separated by ' ' from dictionary.tables 
where libname='WORK' and nobs = 0; 
quit; 

&deletelist ;
proc datasets nolist; 
delete &deletelist; 
quit;




data allele_frq_minDepth5 ;
set allele_frq_minDP5_:;
if mult_altBase_num = . then mult_altBase_num = 0;
run ;

data allele_freq_minDepth10 ;
set allele_frq_minDP10_:;
if mult_altBase_num = . then mult_altBase_num = 0;
run ;

proc sort data = allele_frq_minDepth5 ;
by isolate_pair ;
run;

proc means data = allele_frq_minDepth5  ;
by isolate_pair ;
run;

output out =minDepth5  mean=mean sum=sum nobs=obs ;
run;




