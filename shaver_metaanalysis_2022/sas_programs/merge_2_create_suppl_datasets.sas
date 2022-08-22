libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;
libname manu "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata/manuscript"  ;


/*

create supplemental files 

rename TCA = CM
rename NI = NS

    mass spec:
        for each tech:  (1) merge MA results (strain and path), all annotations and sirius output
                        (2) PH analysis ready
                        (3) export design files with pair info
        
        inputs: work.ma_&tech._&strain 
                work.ma_&tech._&path
                work.ph_&tech.
                work.sr_&tech. 
                dsgn_GT_RP_NEG_pairs_slaw.tsv
                dsgn_GT_RP_POS_pairs_slaw.tsv
                dsgn_GT_HILIC_POS_pairs_slaw.tsv
  nmr
        for each tech:  (1) merge MA results (strain and path)
                        (2) ppm analysis ready
                        (3) export design files with pair info

        inputs: work.ma_&tech._&strain.
                work.ma_&tech._&path.
                work.ph_&tech.
                nmr/dsgn_cdcl3_sbys.tsv
                nmr/dsgn_d2o_sbys.tsv



*/


/* merge by tech --> Mass Spec datasets */
%macro byTech (tech) ;

data &tech._ma ;
merge ma_&tech._: (in=in1) an_&tech._: sr_&tech ph_&tech.;
by featureID ;
if in1 ;
run ;
%mend ;

%byTech (rp_neg);   /* 377 obs */
%byTech (rp_pos);   /*  3953 obs */
%byTech (hilic_pos);/* 199 obs */

/* checking ....
proc contents data = hilic_pos_all ; run ;
proc freq data = hilic_pos_all ;
tables flag_sig_AUM2073 * flag_sig_CB4856 * flag_sig_CX11314 * flag_sig_DL238 * flag_sig_KJ550 * flag_sig_N2 flag_sig_RB2011 * flag_sig_RB2055 * flag_sig_RB2347 * flag_sig_RB2550 * 
flag_sig_UGT49 * flag_sig_UGT60 * flag_sig_VC1265 * flag_sig_VC2524 / out = aaa ;
run ;
*/


/* merge by tech --> NMR datasets */
%macro byTech2 (tech) ;

data &tech._ma ;
merge ma_&tech._: ;
by featureID ;
run ;
%mend ;

%byTech2 (d2o);  /* 585 obs */
%byTech2 (cdcl3);/* 487 obs */



/* save sas datasets and export tsv */
%macro saving (tech, type) ;

data manu.suppl_&tech._meta_table ;
set &tech._ma ;
rename effect_NI = effect_NS ; rename pval_NI = pval_NS ; rename flag_sig_NI = flag_sig_NS;
rename effect_TCA = effect_CM ; rename pval_TCA = pval_CM ; rename flag_sig_TCA = flag_sig_CM ;
run;

proc export data = manu.suppl_&tech._meta_table
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/manuscript/supplemental_data/suppl_&tech._meta_table.tsv"
dbms = tab replace ;
run ;

/* PH/ppm data */
data manu.suppl_&tech._&type._table ;
set ph_&tech. ;
run ;

proc export data = manu.suppl_&tech._&type._table
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/manuscript/supplemental_data/suppl_&tech._&type._table.tsv"
dbms = tab replace ;
run ;

%mend ;

%saving (rp_neg, ph) ;
%saving (rp_pos, ph) ;
%saving (hilic_pos, ph) ;

%saving (cdcl3, ppm) ;
%saving (d2o, ppm) ;


/* export design files as well -- mass spec (all techs are ID!) */
%macro designs (tech, datain, path);

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/&path./&datain."
out = design_&tech. 
dbms = tab replace ;
run;

data manu.design_mass_spec ;
set design_&tech.;
run ;

proc export data = manu.design_mass_spec
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/manuscript/supplemental_data/design_mass_spec.tsv"
dbms = tab replace ;
run ;

%mend ;

%designs (rp_neg, dsgn_GT_rp_neg_pairs_slaw.tsv, design_files) ;


/* design nmr --> cdcl3 and d20 are different!! */
%macro designs (tech, datain, path);

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/&path./&datain."
out = design_&tech. 
dbms = tab replace ;
run;

data manu.design_nmr_&tech. ;
set design_&tech. ;
drop oldSampleID ;
run ;

proc export data = manu.design_nmr_&tech.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/manuscript/supplemental_data/design_nmr_&tech..tsv"
dbms = tab replace ;
run ;

%mend ;


%designs (cdcl3, dsgn_cdcl3_sbys.tsv, nmr) ;
%designs (d2o, dsgn_d2o_sbys.tsv, nmr) ;


/* create list of meta table columns */
proc contents data = manu.suppl_rp_pos_meta_table out = column_name (keep=varnum name); run ;

data names ;
set column_name ;
length varName $52 ;
if find(name, "aos") ge 1 then do ;
    keep = substr(name,1,7);
    varName = compress(keep||"num") ; end ;

if find(name, "effect") ge 1 then do ;
    keep2 = substr(name,1,7) ;
    varName = compress(keep2||"strain/pathway") ; end ;

if find(name, "pval") ge 1 then do ;
    keep3 = substr(name,1,5);
    varName = compress(keep3||"strain/pathway") ; end ;

if find(name, "flag_sig") ge 1 then do ;
    keep4 = substr(name,1,9) ;
    varName = compress(keep4||"strain/pathway") ; end ;

if find(name, "pool_mutant") ge 1 then do ;
    keep5 = substr(name,1,14) ;
    varName = compress(keep5) ; end ;
if find(name, "pool_pd1074") ge 1 then do ;
    keep5 = substr(name,1,14) ;
    varName = compress(keep5) ; end ;
if find(name, "pool_natural") ge 1 then do ;
    keep5 = substr(name,1,15) ;
    varName = compress(keep5) ; end ;

if varName = "" then varName = name ;
keep varName varNum ;
run ;

proc sort data = names ;
by varName ;
run ;

data names2 ;
set names ;
by varName ;
if first.varName then count = 0;
    count +1 ;
if first.varName ;
run ;


proc sort data = names2 ;
by varNum ;
run ;

data meta_table_vars ;
retain varName description ;
length description  $256. ;
set names2 ;

if varName = "featureID" then description = "unique feature identifier" ;
else if find(varName, "effect") ge 1 then description = "meta-analysis effect size for indicated strain or pathway" ;
else if find(varName, "pval") ge 1 then description = "meta-analysis p-value for indicated strain or pathway" ;
else if find(varName, "flag_sig") ge 1 then description = "binary indicator variable. Equal to 1 if p-value <= 0.05" ;

else if varNum ge 90 and varNum le 104 then description = "see Sirius reference for column descriptors" ;
else if find(varName, "aos_num") ge 1 then description = "batch number followed by sample number, see design file for sample information" ;
else if find (varName, "pool") ge 1 then description = "batch number followed by pool type" ;

else description = " " ;
if varName = "rank" then delete ;
drop count  varNum ;
run ;

data manu.meta_table_columns ;
set meta_table_vars ;
run ;

proc export data = manu.meta_table_columns
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/manuscript/supplemental_data/meta_table_columns.tsv"
dbms = tab replace ;
run ;
