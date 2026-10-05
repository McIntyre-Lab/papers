

%macro thr (thr) ;
%macro imp (species) ;

proc import datafile = "!MCLAB/sex_specific_splicing/submission/crossTab_exonERP_vs_sexClass_&species._thr&thr..csv"
out = a_&thr._&species. 
dbms = csv replace ;
guessingrows = MAX ;
run ;

data b_&thr._&species. ;
set a_&thr._&species.;
species = "&species." ;
drop F_bias M_bias no_ttest ;
run;

proc transpose data = b_&thr._&species. out = tall_&thr._&species. ;
by species exonERP ;
run ;

data tall2_&thr._&species. ;
set tall_&thr._&species. ;
label _name_ = "bias_status" ;
rename _name_ = bias_status ;
rename col1 = count ;
run ;

ods output CMH=cmh_&thr._&species.
           ChiSq=chisq_&thr._&species.
           CrossTabFreqs=freq_&thr._&species.;

proc freq data=tall2_&thr._&species.;
    tables exonERP*bias_status / cmh CHISQ;
    weight count;
run;
ods output close ;

/* keep row mean scores from cmh */
data cmh2_&thr._&species. ;
set cmh_&thr._&species.;
where statistic = 2 ;
species = "&species.";
run ;

%mend ;
%imp (dmel6) ;
%imp (dsim2) ;
%imp (dser1) ;
%mend ;
%thr (2_7_12) ;
%thr (2_4_7) ;

data cmh_output_thr2_7_12 ;
retain species ;
set cmh2_2_7_12_: ;
run;

data cmh_output_thr2_4_7 ;
retain species ;
set cmh2_2_4_7_: ;
run;



proc export data = cmh_output_thr2_7_12 
outfile = "!MCLAB/sex_specific_splicing/submission/CMH_output_exonERP_vs_sexClass_thr2_7_12.csv"
dbms = csv replace ;
run;

proc export data = cmh_output_thr2_4_7 
outfile = "!MCLAB/sex_specific_splicing/submission/CMH_output_exonERP_vs_sexClass_thr2_4_7.csv"
dbms = csv replace ;
run;



