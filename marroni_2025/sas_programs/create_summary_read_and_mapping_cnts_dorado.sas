
filename mymacros "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);


/*
QC for 01h_19july2022 run after basecalling with dorado


create table with counts 
  starting read nums  -- raw off machine
  pychopper read nums  *** note pychopper may split reads based on 'fused' reads!!
  mapping read nums (does NOT inlude unmapped and secondary and supplementary)
 
Sum fl_rsc_mapped_read_num + unclass_mapped_read_num


Column headers in output file:
sampleID
start_read_cnts
numReads_GR1kb
perReads_GR1kb
numReads_GR2kb
perReads_GR2kb
ave_RL
med_RL
longestRead
shortestRead
num_fl
num_rescue
num_unclass

fl_rsc_input_2_aln_read_num
fl_rsc_mapped_read_num

unclass_input_2_aln_read_num
unclass_mapped_read_num

sum_map_fl_rsc_unc

*/

/* import design file with sampleID */
/*proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv"
out = df
dbms = csv replace ;
guessingrows = max ;
run ; 
*/

    data WORK.DF    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv" 
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample $16. ;
       informat sampleID $20. ;
       informat sample_w_HL $21. ;
       informat genotype $7. ;
       informat age $3. ;
       informat high_low $4. ;
       informat rep best32. ;
       informat runID best32. ;
       informat TR $3. ;
       informat tech $6. ;
       informat newDate $12. ;
       informat runName $22. ;
       informat ont_dir $113. ;
       informat bcnum best32. ;
       format sample $16. ;
       format sampleID $20. ;
       format sample_w_HL $21. ;
       format genotype $7. ;
       format age $3. ;
       format high_low $4. ;
       format rep best12. ;
       format runID best12. ;
       format TR $3. ;
       format tech $6. ;
       format newDate $12. ;
       format runName $22. ;
       format ont_dir $113. ;
       format bcnum best12. ;
    input
               sample  $
               sampleID  $
               sample_w_HL  $
               genotype  $
               age  $
               high_low  $
               rep
               runID
               TR  $
               tech  $
               newDate $
               runName  $
               ont_dir  $
               bcnum
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

data design ;
retain sample sample_w_HL sampleID TR tech ;
length tech $6. ;
set df ;
TR = compress('TR'||runID) ;
sampleID = compress(sample||'_TR'||runID) ;
*keep sampleID ;
if find(runName, "prom") ge 1 then tech = "prom" ;
else if find(runName, "first") ge 1 then tech = "minion" ;
else if find(ONT_dir, "McIntyre") ge 1 then tech = "prom" ;
else if find(ONT_dir, "_FAR") ge 1 then tech = "minion" ;
else tech = "oops" ;
run ;

/* 01h_19july2022 */
data design_sampleID ;
set design ;
if age = "01h" and newDate = "19july2022" ;
keep sampleID ; 
run;

proc sort data = design_sampleID nodups ;
by _all_ ;
run;

data test ;
set design_sampleID (obs = 3) ;
run;


/* import starting read nums */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/dorado_basecalling/read_cnts_dorado_1to1.csv"
out = start_nums
dbms = csv replace ;
guessingrows = max ;
run ;

data start_nums2 ;
set start_nums ;
rename numReadsGR1kb = numReads_GR1kb ;
rename perReadsGR1kb = perReads_GR1kb ;
run ;

/* import mapped read counts and create single row for each sample */
%macro imp_out (type) ;

%macro imp (sampleID) ;

    data WORK.M_&sampleID._&type.    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/TB14/TB14/pxs_lmm_dros_data/dorado_basecalling/mapped_read_cnts/&sampleID._&type._mapped_read_cnts.csv"
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample $36. ;
       informat input_read_num 32. ;
       informat mapped_read_num 32. ;
       informat prop_mapped 32. ;
       format sample $36. ;
       format input_read_num 32. ;
       format mapped_read_num 32. ;
       format prop_mapped 32. ;
    input
                sample  $
                input_read_num
                mapped_read_num
                prop_mapped
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;
/*
proc import datafile = "/home/ammorse/pxs_lmm_dros_data/mapped_cnts/&sampleID._&type._mapped_read_cnts.csv"
out = m_&sampleID._&type.
dbms = csv replace ;
run ;
*/
data m2_&sampleID._&type. ;
retain sampleID ;
set m_&sampleID._&type.;
sampleID = tranwrd(sample, "_&type", "");
rename input_read_num = &type._input_2_aln_read_num ;
rename mapped_read_num = &type._mapped_read_num ;
drop sample prop_mapped ;
run;

proc sort data = m2_&sampleID._&type. ;
by sampleID ;
run;

%mend ;

%iterdataset(dataset=design_sampleID, function=%nrstr(%imp(&sampleID);));
/*%imp  (dm11037_01h_rep1_TR2); */

%mend ;
%imp_out (fl_rsc) ;
%imp_out (unclass) ;


/* import pychopper counts and create single row for each sample */
%macro imp_out2 (type) ;

%macro imp2 (sampleID) ;

proc import datafile = "/TB14/TB14/pxs_lmm_dros_data/dorado_basecalling/pychop_cnts/&sampleID._&type._pychop_read_cnts.csv"
out = p_&sampleID._&type.
dbms = csv replace ;
guessingrows = max ;
run ;

proc sort data = p_&sampleID._&type. ;
by sampleID ;
run;

%mend ;

%iterdataset(dataset=design_sampleID, function=%nrstr(%imp2(&sampleID);));
/*%imp2  (dm11037_01h_rep1_TR2); */

%mend ;
%imp_out2 (fl) ;
%imp_out2 (rescue) ;
%imp_out2 (unclass) ;



/* create table with starting counts, pychopper counts and mapping counts for each sample */
data table ;
merge start_nums2 p_: m2_: ;
by sampleID ;
run;

proc contents data  = table ; run;

data table2 ;
retain sampleID start_read_cnts numReads_GR1kb perReads_GR1kb numReads_GR2kb perReads_GR2kb ave_RL med_RL longestRead shortestRead num_fl num_rescue num_unclass 
    sum_pychop_fl_rsc_unclass sum_pychop_fl_rsc ;
set table ;
where num_fl ne . ;
sum_pychop_fl_rsc_unclass = (num_fl + num_rescue + num_unclass) ;
sum_pychop_fl_rsc = (num_fl + num_fl + num_rescue) ;

sum_map_fl_rsc_unc = (fl_rsc_mapped_read_num + unclass_mapped_read_num) ;

run;

/* merge in design info */
data dsgn ;
retain sample sampleID genotype age rep high_low sample_w_HL TR tech newDate ;
set design ;
keep sample sampleID genotype rep age high_low sample_w_HL TR tech newDate  ;
run ;

proc sort data = dsng ;
by sampleID ;
proc sort data = table2 ;
by sampleID ;
run ;

data table_w_dsgn ;
merge dsgn (in=in1) table2 (in=in2) ;
by sampleID ;
if in2 ;
run;


proc export data = table_w_dsgn 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/patrika_rils/dorado_basecalling/read_cnt_and_mapping_cnt_table_qc_dorado.csv"
dbms = csv replace ;
run ;




