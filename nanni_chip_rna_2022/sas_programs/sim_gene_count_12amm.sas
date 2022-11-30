
libname chiprna '!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs';


data sim_chip_rna_frag_flags_anno;
set   chiprna.sim_chip_rna_frag_flags_anno;
run;
;

proc contents data = sim_chip_rna_frag_flags_anno ; run;


proc sort data=sim_chip_rna_frag_flags_anno ;
by fbgn;
run;

proc sort data=count_transcripts1;
by fbgn;
run;

/* dump so don't accidentally grab mel files */
proc datasets library = WORK ;
delete flip_: 
count_2ratio
count_expression
count_k4
count_k27
count_rna
count_rna05
gene_:
;
run;


/*at the gene level eliminate teh multigene*/
proc freq data=sim_chip_rna_frag_flags_anno noprint;
by fbgn;
where featuretype="fragment" and flag_multigene=0 ;
tables sex_bias/out=count_2ratio;
tables ratio_expressed/out=count_ratio2;
tables ratio_trend/out=count_trend; *for trend in expression;
tables k4_detected/out=count_k4;
tables k27_detected/out=count_k27;
tables rna_detected/out=count_rna;
tables flag_rna_expressed/out=count_rna05; *new;
run;


/* new piece expressed at apn05*/

proc transpose data=count_rna05 out=flip_rna05;
by fbgn;
var count;
id flag_rna_expressed;
run;

data gene_rna05;
set flip_rna05;

if _1=. then _1=0;
if _0=. then _0=0;
rename _1=num_frag_detected05
		_0=num_frag_not_detected05;
sum_fragments_rna_detected05=_0+_1;
if _0<sum_fragments_rna_detected05 then flag_rna_detected05=1;
	else flag_rna_detected05=0;
drop _name_ _label_;
	run;

proc freq data=gene_rna05;
tables flag_rna_detected05;
run;


/* new piece ratio_fragments detected and at 2x*/
/*note that the original dat set contains missing labeles indicating fragments present but not detected
    6128 obs omitted due to missing ID values (no ratio_trend)  */

proc transpose data=count_trend out=flip_trend;
by fbgn;
var count;
id ratio_trend;
run;

data check ;
retain fbgn ratio_trend ;
set sim_chip_rna_frag_flags_anno;
where fbgn = "FBgn0012820";
run;  /* these are all f-m ambiguous */



/*6128 values ommitted due to missing values this should be ok*/
/*perhaps check in previous program?*/

data gene_ratio_trend;
set flip_trend;

rename male=num_ratio_male fem=num_ratio_fem ;

if male=. then male=0;
if fem=. then fem=0;

sum_fragments_ratio_trend=male + fem;

if male > 0 and fem > 0 
			then gene_ratio_trend="male_and_female";
	else if male>0 
			then gene_ratio_trend="male";
	else if fem>0  
			then gene_ratio_trend="fem";
	else gene_ratio_trend="oops";

if sum_fragments_ratio_trend=0 then gene_ratio_trend="none";
drop _name_ _label_;
	run;

proc freq data=gene_ratio_trend;
tables gene_ratio_trend;
run; /*
gene_ratio_
trend              Frequency
----------------------------
fem                    2937
male                   4208
male_and_female        5121
none                   2880  */

proc transpose data=count_ratio2 out=flip_ratio2;
by fbgn;
var count;
id ratio_expressed;
run;

/*6128 values ommitted due to missing values this should be ok*/
/*   AMM looked at in prev program, these are missing or female-male bias not consistent, or multigene
     looks ok*/

data gene_ratio2;
set flip_ratio2;

rename male=num_ratio2_male fem=num_ratio2_fem ;

if male=. then male=0;
if fem=. then fem=0;
if unb=. then unb=0 ;

sum_fragments_ratio2=male+fem+unb;

if male>0 and fem>0 
			then gene_ratio2="male_and_female";
	else if male>0 
			then gene_ratio2="male";
	else if fem>0  
			then gene_ratio2="fem";
    else if unb>0 
            then gene_ratio2="";  /* make unbiased frag count missing */
	else gene_ratio2="oops";

if sum_fragments_ratio2=0 then gene_ratio2="none";
drop _name_ _label_;
	run;

proc freq data=gene_ratio2;
tables gene_ratio2;
run; /*
gene_ratio2        Frequency
----------------------------
fem                     463
male                   1746
male_and_female          35
none                   2880

freq missing = 10022  */


/* end new code except for check piece after next segment*/

/*detected*/
proc transpose data=count_rna out=flip_rna;
by fbgn;
var count;
id rna_detected;
run;

data gene_rna;
set flip_rna;

rename both=num_rna_both
male=num_rna_male
fem=num_rna_fem 
none=num_rna_none;

if both=. then both=0;
if male=. then male=0;
if fem=. then fem=0;
if none=. then none=0;

sum_fragments_rna_detected=male+fem+both;

if male>0 and fem>0  and both>0
			then gene_rna="male_and_female_and_both";
	else if male>0 and fem>0  
			then gene_rna="male_and_female";
	else if male>0 and both>0
			then gene_rna="male_and_both";
	else if fem>0  and both>0
			then gene_rna="fem_and_both";
	else if both>0
			then gene_rna="both";
	else if male>0 
			then gene_rna="male";	
	else if fem>0 
			then gene_rna="fem";
	else gene_rna="oops";

if sum_fragments_rna_detected=0 then gene_rna="none";
drop _name_ _label_;
	run;

proc freq data=gene_rna;
tables gene_rna;
run;


/* check that nothing detected at 05 is missing at apn0 */

proc sort data=gene_rna;
by fbgn;

proc sort data=gene_rna05;
by fbgn;

data check;
merge gene_rna gene_rna05;
by fbgn;
run;

proc freq data=check;
tables flag_rna_detected05*gene_rna;
run;

/*
flag_rna_detected05     gene_rna

Frequency|
Percent  |
Row Pct  |
Col Pct  |both    |fem     |fem_and_|male    |male_and|male_and|male_and|none    |  Total
         |        |        |both    |        |_both   |_female |_female_|        |
         |        |        |        |        |        |        |and_both|        |
---------+--------+--------+--------+--------+--------+--------+--------+--------+
       0 |   2339 |    195 |    332 |    351 |    390 |     27 |    123 |   2826 |   6583
         |  15.44 |   1.29 |   2.19 |   2.32 |   2.57 |   0.18 |   0.81 |  18.66 |  43.46
         |  35.53 |   2.96 |   5.04 |   5.33 |   5.92 |   0.41 |   1.87 |  42.93 |
         |  23.25 | 100.00 |  45.42 |  99.72 |  51.52 | 100.00 |  62.76 | 100.00 |
---------+--------+--------+--------+--------+--------+--------+--------+--------+
       1 |   7723 |      0 |    399 |      1 |    367 |      0 |     73 |      0 |   8563
         |  50.99 |   0.00 |   2.63 |   0.01 |   2.42 |   0.00 |   0.48 |   0.00 |  56.54
         |  90.19 |   0.00 |   4.66 |   0.01 |   4.29 |   0.00 |   0.85 |   0.00 |
         |  76.75 |   0.00 |  54.58 |   0.28 |  48.48 |   0.00 |  37.24 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+--------+--------+
Total       10062      195      731      352      757       27      196     2826    15146
            66.43     1.29     4.83     2.32     5.00     0.18     1.29    18.66   100.00


/*k4*/
proc transpose data=count_k4 out=flip_k4;
by fbgn;
where k4_detected ne "";
var count;
id k4_detected;
run;

data gene_k4;
set flip_k4;


rename both=num_k4_both
male=num_k4_male
fem=num_k4_fem 
none=num_k4_none;

if both=. then both=0;
if male=. then male=0;
if fem=. then fem=0;
if none=. then none=0;

sum_fragments_k4=male+fem+both;

if male>0 and fem>0  and both>0
			then gene_k4="male_and_female_and_both";
	else if male>0 and fem>0  
			then gene_k4="male_and_female";
	else if male>0 and both>0
			then gene_k4="male_and_unb";
	else if fem>0  and both>0
			then gene_k4="fem_and_unb";
	else if both>0
			then gene_k4="unb";
	else if male>0 
			then gene_k4="male";	
	else if fem>0 
			then gene_k4="fem";
	else gene_k4="oops";



if male>0 and fem>0  and both>0
			then gene_k4_sex="male_and_female";
	else if male>0 and fem>0  
			then gene_k4_sex="male_and_female";
	else if male>0 and both>0
			then gene_k4_sex="male";
	else if fem>0  and both>0
			then gene_k4_sex="fem";
	else if both>0
			then gene_k4_sex="unb";
	else if male>0 
			then gene_k4_sex="male";	
	else if fem>0 
			then gene_k4_sex="fem";
	else gene_k4_sex="oops";

if gene_k4_sex="male" or gene_k4_sex="fem" or gene_k4_sex="male_and_female" then gene_k4_bias="sex_bias";
else gene_k4_bias=gene_k4_sex;


if sum_fragments_k4=0 then gene_k4="none";
if sum_fragments_k4=0 then gene_k4_sex="none";
if sum_fragments_k4=0 then gene_k4_bias="none";

drop _name_ _label_;
	run;

proc freq data=gene_k4;
tables gene_k4*gene_k4_sex;
tables gene_k4_sex*gene_k4_bias;
run;


/*k27*/
proc transpose data=count_k27 out=flip_k27;
by fbgn;
where k27_detected ne "";
var count;
id k27_detected;
run;

data gene_k27;
set flip_k27;

rename both=num_k27_both
male=num_k27_male
fem=num_k27_fem 
none=num_k27_none;

if both=. then both=0;
if male=. then male=0;
if fem=. then fem=0;
if none=. then none=0;

sum_fragments_k27=male+fem+both;

if male>0 and fem>0  and both>0
			then gene_k27="male_and_female_and_both";
	else if male>0 and fem>0  
			then gene_k27="male_and_female";
	else if male>0 and both>0
			then gene_k27="male_and_unb";
	else if fem>0  and both>0
			then gene_k27="fem_and_unb";
	else if both>0
			then gene_k27="both";
	else if male>0 
			then gene_k27="male";	
	else if fem>0 
			then gene_k27="fem";
	else gene_k27="oops";


if male>0 and fem>0  and both>0
			then gene_k27_sex="male_and_female";
	else if male>0 and fem>0  
			then gene_k27_sex="male_and_female";
	else if male>0 and both>0
			then gene_k27_sex="male";
	else if fem>0  and both>0
			then gene_k27_sex="fem";
	else if both>0
			then gene_k27_sex="unb";
	else if male>0 
			then gene_k27_sex="male";	
	else if fem>0 
			then gene_k27_sex="fem";
	else gene_k27_sex="oops";

if gene_k27_sex="male" or gene_k27_sex="fem" or gene_k27_sex="male_and_female" then gene_k27_bias="sex_bias";
else gene_k27_bias=gene_k27_sex;


if sum_fragments_k27=0 then gene_k27="none";
if sum_fragments_k27=0 then gene_k27_sex="none";
if sum_fragments_k27=0 then gene_k27_bias="none";


drop _name_ _label_;
	run;

proc freq data=gene_k27;
tables gene_k27*gene_k27_sex;
tables gene_k27_sex*gene_k27_bias;
run;

/*combined k4 and rna*/

proc transpose data=count_2ratio out=flip_s;
by fbgn;
where sex_bias ne "";
var count;
id sex_bias;
run;

data classify_mixed;
set flip_s;

if fem_full=. then fem_full=0;
if fem_exp=. then fem_exp=0;
if male_full=. then male_full=0;
if male_exp=. then male_exp=0;
if male_swch=. then male_swch=0;
if fem_swch=. then fem_swch=0;
if unb=. then unb=0;

sum_male_express=male_full+male_exp+male_swch;
sum_female_express=fem_full+fem_exp+fem_swch;
num_fragments_express=fem_full+fem_exp+male_full+male_exp+male_swch
	+fem_swch+unb;

if sum_male_express>0 and sum_female_express>0 then gene_sex_bias="male_and_female";
	else if sum_male_express>0 then gene_sex_bias="male";
	else if sum_female_express>0  then gene_sex_bias="female";
	else gene_sex_bias="unbiased";

if num_fragments_express=0 then gene_sex_bias="";

if fem_full>0 then flag_fem_full=1;
	else flag_fem_full=0;
if fem_exp>0 then flag_fem_exp=1;
	else flag_fem_exp=0;
if male_full>0 then flag_male_full=1;
	else flag_male_full=0;
if male_exp>0 then flag_male_exp=1;
	else flag_male_exp=0;
if male_swch>0 then flag_male_swch=1;
	else flag_male_swch=0;
if fem_swch>0 then flag_fem_swch=1;
	else flag_fem_swch=0;
if unb>0 then flag_unb_k4e=1;
	else flag_unb_k4e=0;



drop _name_ _label_;

	run;


proc freq data=classify_mixed;
tables gene_sex_bias;
run;

proc sort data=classify_mixed;
by fbgn;

proc sort data=gene_rna;
by fbgn;


proc sort data=gene_k4;
by fbgn;


proc sort data=gene_k27;
by fbgn;
run;

data sim_gene_flags_anno;
merge gene_k4 (in=in1) gene_k27 (in=in2) gene_rna (in=in4) gene_rna05 (in=in7) gene_ratio2 gene_ratio_trend classify_mixed (in=in3)
fbgn_anno_all(in=in5) sim_fbgn_2_coord(in=in6) count_transcripts1;
by fbgn;
if in1 or in3;

if chrom="Scf_X" then xsome="X";
	else if chrom="Scf_4" then xsome="4";
	else if chrom="Scf_Y" then xsome="Y";
	else if chrom="Scf_2L" or chrom="Scf_2R" or chrom="Scf_3L" or chrom="Scf_3R" then xsome="A";
	else xsome="";

if chrom="Scf_X" then XCHR=1; 
	else if  xsome="" then xchr="";
	else XCHR=0;
if chrom="Scf_4" or chrom="Scf_Y" then flag_4Y=1;
	else flag_4Y=0;
if chrom="Scf_X" then flag_X=1;
	else flag_X=0;
if chrom="Scf_Y" then flag_Y=1;
	else flag_Y=0;
run;

*fbgn2coord;

proc freq data=sim_gene_flags_anno;
tables chrom;
run;


data chiprna.sim_gene_flags_anno ;
set sim_gene_flags_anno ;
run;


PROC EXPORT DATA= chiprna.sim_gene_flags_anno 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/data_files/sim_gene_flags_anno.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



title "sim gene level single transcripts only";
proc freq data= chiprna.sim_gene_flags_anno;
where flag_single_transcript=1;
tables flag_fem_full*flag_male_full*gene_k4 ;
tables gene_sex_bias*gene_rna;
tables gene_k4*gene_k27;
tables gene_sex_bias gene_rna gene_k4 gene_k27;
tables xsome*(gene_k4 gene_k27 gene_rna gene_sex_bias);
tables xsome*(flag_male_full flag_fem_full flag_male_exp flag_fem_exp
flag_male_swch flag_fem_swch flag_unb_k4e);
run;


title "sim gene level multi transcripts only";
proc freq data= chiprna.sim_gene_flags_anno;
where flag_single_transcript=0;
tables flag_fem_full*flag_male_full*gene_k4 ;
tables gene_sex_bias*gene_rna;
tables gene_k4*gene_k27;
tables gene_sex_bias gene_rna gene_k4 gene_k27;
tables xsome*(gene_k4 gene_k27 gene_rna gene_sex_bias);
tables xsome*(flag_male_full flag_fem_full flag_male_exp flag_fem_exp
flag_male_swch flag_fem_swch flag_unb_k4e);
run;


ods pdf file="!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/sim_gene_results_v2.pdf";


title "sim gene level ";
proc freq data= chiprna.sim_gene_flags_anno;
tables flag_fem_full*flag_male_full*gene_k4 ;
tables gene_sex_bias*gene_rna;
tables gene_k4*gene_k27;
tables gene_sex_bias gene_rna gene_k4 gene_k27 gene_k4_sex gene_k27_sex;
tables xsome*(gene_k4 gene_k27 gene_rna gene_sex_bias gene_k4_sex gene_k27_sex);
tables xsome*(flag_male_full flag_fem_full flag_male_exp flag_fem_exp
flag_male_swch flag_fem_swch flag_unb_k4e);
run;

ods pdf close;


title "sim gene level single vs multi ";
proc freq data= chiprna.sim_gene_flags_anno;
tables flag_single_transcript*gene_sex_bias;
tables flag_single_transcript*xsome*(gene_k4 gene_k27 gene_rna gene_sex_bias);
run;


data FBgn0086038;
set chiprna.sim_chip_rna_frag_flags_anno;;
where fbgn="FBgn0086038";
run;

