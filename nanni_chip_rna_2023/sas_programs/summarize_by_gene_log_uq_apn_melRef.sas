
libname ortho "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq_ortho";
libname drosRNA "!MCLAB/Dros_CHIP_RNA_ms/sas_data/RNAseq";
libname df  "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/svn_lmm_dros_head_data/ethanol_srna/sas_data";


/*

starting with normalized coverage counts on fusions from adalena for sim to mel and mel to mel
    log_uq_apn

     mel anno to cnt on fusions
     roll to gene

APN_gene = sum of all fusions * fusionLength divided by the sum of all fusionLengths
                       = (length1*apn1 + length2*apn2) / (length1+length2)

proc means use weight statement
	weight log_uq_apn
	var length

** rolled to gene level including fusions with 0 reads (weighted log_up_apn) AND omitting fusions with 0 reads (wieghted log_uq_apn_noZero)

Fusion coverage counts:       /blue/mcintyre/share/etoh_srna/rnaseq/coverage_cnts/fusion_coverage_cnts/*
        mel to mel genome alignments are in the format cvrg_cnts_mel_${GENO}_${SEX}_${TRT}_rep${REP}_combined.csv
        sim to mel genome alignments are in the format cvrg_cnts_sim_${GENO}_${SEX}_${TRT}_rep${REP}_sim2mel.csv

sampleID example:  mel_301_m_noEtoh_rep2

Fusion annotation (w/gene_id):  /blue/mcintyre/share/references/dmel_fb617/dmel617_fusion_annotations.csv 
    this is ID to /McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_fusion_annotations.csv

design for samples to loop over:  df.design_ethanol_srna_rnaseq
*/


filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);


/* fusions 2 genes anno  */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_fusion_annotations.csv"
out = fusion2gene 
dbms = csv replace ;
run ;  /* 64619 obs */

data fusion2gene2 ;
set fusion2gene ;
rename fusion_id = feature_ID ;
where flag_multigene = 0 ;
keep fusion_id gene_id ;
run;
 /* 60815 features NOT multigene  */

proc sort data = fusion2gene2 dupout = chk2 nodups ;
by _all_ ;
run;  /* 60,815 (0 dups) */


/* get fusion lengths from bed file */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_fusions_coverage.bed"
out = fusionBed 
dbms = tab replace ;
guessingrows = MAX ;
getnames = no ;
run ;

data fusionLen ;
set fusionBed ;
fusionLength = (var3 - var2 + 1) ;
keep var4 fusionLength ;
rename var4 = feature_id ;
run ;

/* add fusion lengths to fusion2gene */
proc sort data = fusion2gene2 ;
by feature_id ;
proc sort data = fusionLen ;
by feature_id ;
run;

data fusion2gene_len ;
merge fusion2gene2 (in=in1)  fusionLen (in=in2) ;
by feature_id ;
if in1 and in2 ;
run ;  /* 60815 obs */


/* list of samples to loop over format - mel_301_m_noEtoh_rep2 */
data dsgn2 ;
set df.design_ethanol_srna_rnaseq ;
keep species genotype sex arbitrary_rep treatment;
run ;

proc sort data = dsgn2 nodups ;
by _all_ ;
run;

data dsgn ;
retain sampleID ;
set dsgn2 ;
sampleID = compress(species||'_'||genotype||'_'||sex||'_'||treatment||'_rep'||arbitrary_rep) ;
keep sampleID ;
if sampleID = "sim_12_m_noEtoh_rep2"  then delete ;
run;  /* 47 samples */


proc sort data = dsgn nodups ;
by _all_ ;
run; /* uniq */

data test_list ;
set dsgn ;
where sampleID = "sim_11_f_etoh_rep1" ;
run ;


data normed;
set ortho.norm_stats_fusion_melRef ;
rename featureID = feature_id ;
run ; /* 2,362,079 obs */

proc sort data = fusion2gene_len ;
by feature_id ;
proc sort data =normed ;
by feature_id ;
run;


data normed2 ;
merge normed (in=in1) fusion2gene_len (in=in2) ;
by feature_id ;
if in1 ;
if gene_id = "" then delete ;
run;  /* 2,205,287 obs */

proc sort data = normed2 ;
by gene_id feature_id ;
run;




%macro summing (sampleID) ;

data in2_&sampleID.  ;
set normed2 ;
where sampleID = "&sampleID" ;
sampleID = "&sampleID";
run;

/* set aside columns want to add back in after collapsing to gene level - these remain the same */
data anno_&sampleID. ;
set in2_&sampleID. ;
keep gene_id mapped_reads ;
rename gene_id = feature_id ;
run ;

proc sort data = anno_&sampleID. nodups ;
by _all_ ;
run ;

proc sort data =anno_&sampleID.;
by feature_id ;
run ;

proc sort data = in2_&sampleID. ;
by gene_id ;
run;

/* output in stacked format - for depth and length only */
proc means data = in2_&sampleID. sum stackods ;
by gene_id;
ods output summary = want_&sampleID. ;
run ;

proc means data = in2_&sampleID. sum stackods ;
by gene_id;
where apn ne 0 ;
ods output summary = want_no0_&sampleID. ;
run ;


/* 
/* get list of variable names containing APN values
	for looping through proc sql to calc weighed APN */
proc contents data = in2_&sampleID out = var_&sampleID ; run;

data var2_&sampleID. ;
set var_&sampleID ;
if find(name, "log_uq_apn") ge 1 ;
num = _n_ ;
keep name num ;
run;


/* calc weighted ave for APN values 
	(length1*apn1 + length2*apn2) / (length1+length2)

    *** note this weighted APN SHOULD include any fusions with 0 counts
*/
%macro weighted (name, num) ;

proc sql ;
create table wt_&sampleID._&num. as 
	select 
		gene_id,
                sum(&name.*fusionLength) / sum(fusionLength) as &name 
	from ( 
		select
		*
		from in2_&sampleID
		group by gene_id)
	group by gene_id;
	select * from wt_&sampleID._&num.;
quit ;

data w_&sampleID._&num.;
set wt_&sampleID._&num. ;
keep gene_id &name ;
rename gene_id = feature_ID ;
run ;

proc sort data = w_&sampleID._&num. ;
by feature_ID ;
run ;

%mend ;

%iterdataset(dataset=var2_&sampleID, function=%nrstr(%weighted(&name, &num);));
/* %iterdataset(dataset=var2_mel_W55_M_rep3, function=%nrstr(%weighted(&name, &num);)); */


/* calc weighted ave for APN values --- OMIT FUSIONS WITH 0 COUNTS
	(length1*apn1 + length2*apn2) / (length1+length2)

    *** note this weighted APN should NOT include any fusions with 0 counts
*/


%macro weighted_no0 (name, num) ;

data in3_&sampleID ;
set in2_&sampleID. ;
if log_uq_apn = . then delete ;
run;

proc sql ;
create table wt_no0_&sampleID._&num. as 
	select 
		gene_id,
                sum(&name.*fusionLength) / sum(fusionLength) as &name 
	from ( 
		select
		*
		from in3_&sampleID
		group by gene_id)
	group by gene_id;
	select * from wt_no0_&sampleID._&num.;
quit ;

data w_no0_&sampleID._&num.;
set wt_no0_&sampleID._&num. ;
rename gene_id = feature_ID ;
keep gene_id &name ;
rename &name = &name._noZero ;
run ;

proc sort data = w_no0_&sampleID._&num. ;
by feature_ID ;
run ;

%mend ;

%iterdataset(dataset=var2_&sampleID, function=%nrstr(%weighted_no0(&name, &num);));


/* prep the count data  */
data x2_&sampleID. ;
set want_&sampleID. ;
where variable = "region_depth"  or variable = "region_length" ;
rename gene_id = feature_ID ;
run ;

proc transpose data = x2_&sampleID out = cnt_&sampleID  ;
by feature_id ;
var sum ;
id variable ;
run ;

proc sort data = cnt_&sampleID. ;
by feature_id ;
run ;

data x2_no0_&sampleID. ;
set want_no0_&sampleID. ;
where variable = "region_depth"  or variable = "region_length" ;
rename gene_id = feature_ID ;
run ;

proc transpose data = x2_no0_&sampleID out = cnt_no0_&sampleID  ;
by feature_id ;
var sum ;
id variable ;
run ;

data cnt_no0_&sampleID. ;
set cnt_no0_&sampleID. ;
rename region_depth = region_depth_noZero ;
rename region_length = region_length_noZero ;
run;

proc sort data = cnt_no0_&sampleID. ;
by feature_id ;
run ;


/* merge count data and weighted APN values (already sorted) */
data cc_gene_&sampleID.  ;
merge anno_&sampleID.  cnt_&sampleID. cnt_no0_&sampleID. w_&sampleID._: w_no0_&sampleID._: ;
by feature_id ;
run ;

data cc2_gene_&sampleID. ;
retain feature_id mapped_reads region_depth region_length log_uq_apn region_length_noZero region_depth_noZero ;
format sampleID $26. ;
format region_depth 26. ;
format region_depth_noZero 26. ;
format region_length 26. ;
format region_length_noZero  26. ;  
set cc_gene_&sampleID. ;
rename region_depth = region_depth_mpile ;
sampleID = "&sampleID";
drop _name_ region_depth_noZero ;
run;

data cvrg_gene_&sampleID ;
retain sampleID feature_id ;
set cc2_gene_&sampleID ;
run;

%mend ;

/*%iterdataset(dataset=test_list, function=%nrstr(%summing(&sampleID);)); */
%iterdataset(dataset=dsgn, function=%nrstr(%summing(&sampleID);));  



/* set samples together and save */
data ortho.cvrg_gene_wt_log_uq_apn_melRef ;
set cvrg_gene_:;
rename feature_id = geneID ;
run ;

proc freq data = ortho.cvrg_gene_wt_log_uq_apn_melRef ;
tables sampleID  / out = ck_sampleID ;
run;
    /* 47 sampleIDs, 11,067 obs each */


