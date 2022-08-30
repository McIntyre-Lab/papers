libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

create table containing effect sizes (by strain) for all strains (col) and all features (add suffix to featureID for each tech) include nmr and mass spec



*/


filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



data list ;
set sweet16.dsgn_gt_rp_neg_slaw ;
keep strain ;
if strain = "" then delete ;
if strain = "PD1074" then delete ;
run ;

proc sort data = list  nodups  ;
by _all_ ;
run;

/* import mass spec effect sizes */
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./MA_FE_rank_byMut_&tech._&strain._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data a_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
newID = compress("&tech"||'_'||featureID) ;
drop featureID ;
rename newID = featureID ;
keep effect newID ;
run ;

proc sort data =a_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=list, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (rp_neg) ;
%imp_outer (rp_pos) ;
%imp_outer (hilic_pos) ;


/* import NMR effect sizes */
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_fixedEffect_NMR_&tech._rank_&strain._&tech._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data a_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
newID = compress("&tech"||'_'||featureID) ;
drop featureID ;
rename newID = featureID ;
keep effect newID ;
run ;

proc sort data =a_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=list, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (d2o) ;
%imp_outer (cdcl3) ;


/* merge all strains for each tech */

%macro byTech (tech) ;

data &tech._all ;
merge a_&tech._: ;
by featureID ;
run ;
%mend ;

%byTech (rp_neg);
%byTech (rp_pos);
%byTech (hilic_pos);

%byTech (d2o);
%byTech (cdcl3);

data big_one ;
retain featureID ;
length featureID $26.;
set rp_neg_all rp_pos_all hilic_pos_all d2o_all cdcl3_all ;
run ;


proc export data = big_one 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/analysis_massSpec_and_NMR/effectSizes_MassSpec_NMR.tsv"
dbms = tab replace ;
run;

/* mmc needs design */
proc contents data = big_one out = variables noprint;
run;

data dsgn ;
set variables ;
rename name = sampleID ;
keep name ;
if find(name, "featureID") ge 1 then delete ;
run;

proc export data = dsgn 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/analysis_massSpec_and_NMR/dsgn_effectSizes_MassSpec_NMR.tsv"
dbms = tab replace ;
run;



/* flip so columns are metabolites and rows are strains */
proc transpose data = big_one out = flip NAME =strain_effectSize ;
id featureID ;
run;

data flip2 ;
set flip ;
label strain_effectSize = "strain_effectSize" ;
run ;

proc sort data = flip2 ;
by strain_effectSize ;
run ;

/* export with > 5k variables not working.... 
proc export data = flip2 (obs=0)
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/analysis_massSpec_and_NMR/effectSizes_flipped_MassSpec_NMR.tsv"
dbms=tab replace ;
putnames = yes ;
run ;

data _null_ ;
set flip2 ;
file "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/analysis_massSpec_and_NMR/effectSizes_flipped_MassSpec_NMR.tsv"
mod dlm=tab lrecl=2000000;
put (_all_) (:) ;
run ;
*/

/* mmc needs design */
proc contents data = flip2 out = flip_vars noprint;
run;

data dsgn2 ;
set flip_vars ;
rename name = sampleID ;
keep name ;
if find(name, "strain") ge 1 then delete ;
run;

proc export data = dsgn2 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/analysis_massSpec_and_NMR/dsgn_effectSizes_flipped_MassSpec_NMR.tsv"
dbms = tab replace ;
run;


