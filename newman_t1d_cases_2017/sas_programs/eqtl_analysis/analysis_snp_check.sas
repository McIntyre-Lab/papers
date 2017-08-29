libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* NEED TO CHECK:
1. SNPs in Onengut Supp Table 1 are in LD or tagged SNP lists
2. Kept SNPs are in LD/tagged SNP lists
3. Kept SNPs are in the list of tested SNPs
Then we need to merge these so that we can relate T1D-assoc SNP to a tested eQTL SNP
*/
/* Flag tested SNPs as being T1D gene or not. This is a check to see if at least all T1D SNPs are in the other lists */

data flag_diabetes_snps;
   set eqtl.eqtl_results_summary_table;
   keep snp_id flag_diabetes_gene;
   rename flag_diabetes_gene=flag_diabetes_snp;
run;

proc sort data=flag_diabetes_snps nodup;
   by snp_id flag_diabetes_snp;
run;

*check to see if SNPs are unique -- if not, only take diabetes-flagged duplicate;
proc freq data=flag_diabetes_snps noprint;
   tables snp_id / out=snp_count;
run;

data dup_snps;
   set snp_count;
   if count gt 1;
run;

/* Resort SNP lists, take only the diabetes-flagged duplicate */

proc sort data=flag_diabetes_snps;
   by snp_id descending flag_diabetes_snp;
run;

data tested_snp_list;
   set flag_diabetes_snps;
   by snp_id;
   if first.snp_id then output;
run;

/* Make this permenant in case we want to look at it again later */

data eqtl.tested_snp_list;
   set tested_snp_list;
run;


/* Map tested SNP to LD SNPs */

data ld_snps;
   set eqtl.snp_ld_data;
   rename snp_id_1=snp_id;
run;

proc sort data=ld_snps;
   by snp_id;
proc sort data=eqtl.tested_snp_list;
   by snp_id;
run;

data tested_snps_w_ld;
   merge eqtl.tested_snp_list (in=in1) ld_snps (in=in2);
   by snp_id;
   if in1 and in2 then do;
      flag_tested_snp=1; output; end;
   else if in1 then do;
      snp_id_2=snp_id; r2=1; flag_tested_snp=1; output; end;
   else do;
      flag_tested_snp=0; output; end;
run;

/* Map tested SNP to tagging SNPs */

data tagging_snps;
   set eqtl.tagged_snp_list;
   if tagged_snps='NONE' then tagged_snps=snp_id; *change the "none" to the tagging SNP itself;
run;

proc sort data=tagging_snps;
   by snp_id;
proc sort data=tested_snps_w_ld;
   by snp_id;
run;

data tested_snps_w_tags no_tags no_tests;
   merge tested_snps_w_ld (in=in1) tagging_snps (in=in2);
   by snp_id;
   if in1 and in2 then output tested_snps_w_tags;
   else if in1 then output no_tags; *0 0bs!;
   else output no_tests;  *0 0bs!;
run;
* all tested SNPs have tags;

/* Now I need to stack (un-concatenate) the tagging snps */
data tested_snps_w_tags_2;
    length tagged_snp_id $15.; 
    set tested_snps_w_tags; 
    do i=1 by 1 while(scan(tagged_snps,i,'|') ^=' ');
        tagged_snp_id=scan(tagged_snps,i,'|'); 
        drop i tagged_snps;
        output; 
        end; 
run;

* append the tagging SNP to the list of tagged SNPs;
data tested_snps_w_tags_3;
   set tested_snps_w_tags_2;
   snp_id2=snp_id;
   drop tagged_snp_id;
   rename snp_id2=tagged_snp_id;
run;

proc sort data=tested_snps_w_tags_3 nodup;
   by snp_id tagged_snp_id snp_id_2;
run;

* stack datasets ;
data tested_snps_w_tags_4;
  set tested_snps_w_tags_3 tested_snps_w_tags_2;
run;

* remove duplicates ;
proc sort data=tested_snps_w_tags_4 nodup;
   by snp_id tagged_snp_id snp_id_2;
run;

/* Now only need the tagging SNP data, will drop extraneous LD data and remove duplicate entries */

data tested_snps_w_tags_5;
   set tested_snps_w_tags_4;
   keep flag_diabetes_snp snp_id flag_tested_snp tagged_snp_id;
run;

proc sort data=tested_snps_w_tags_5 nodup;
  by snp_id tagged_snp_id;
run;


/* rename variable and merge in supp table 1 SNPs */
data tested_snps_w_tags_6;
   set tested_snps_w_tags_5;
   rename snp_id=tested_snp_id
          tagged_snp_id=snp_id;
run;

* adding the UCSC-formatted SNP ids for later checking ;
data onengut_snp_list;
  length ucsc_snp_id $15.;
  length snp_id $15.; format snp_id $15.; informat snp_id $15.;
  set eqtl.supptable1_snp_list;
  ucsc_snp_id=catx(":",chr,position);
  keep ucsc_snp_id snp_id;
  run;

proc sort data=onengut_snp_list nodup;
   by snp_id  ucsc_snp_id;
proc sort data=tested_snps_w_tags_6;
   by snp_id;
run;

data onengut_snp_list_tested no_test no_onengut;
   merge onengut_snp_list (in=in1) tested_snps_w_tags_6 (in=in2);
   by snp_id;
   if in1 and in2 then output onengut_snp_list_tested;
   else if in1 then output no_test;
   else output no_onengut;
run;

*2027 SNPs in the Onengut list, 160763 SNPs in the tested-to-tagged list;
* 1768 Onengut SNPs have a corresponding test; * 896 Onengut SNPs do not have a corresponding test. This is okay!;
* 158995 tested-to-tagged SNPs do not have an Onengut list counterpart. This is okay!;

/* Need to also check if UCSC-style "chrZ:WWWWWWWWW" formatted SNPs aren't being dropped */
data ucsc_snps_to_check;
   set no_test;
   keep snp_id ucsc_snp_id; 
   rename snp_id=rs_snp_id
          ucsc_snp_id=snp_id;
run;

data remaining_tests;
   set no_onengut;
   drop ucsc_snp_id;
run;
 
proc sort data=ucsc_snps_to_check;
    by snp_id;
proc sort data=remaining_tests;
    by snp_id;
run;

data ucsc_snps_tested no_test no_ucsc;
   merge ucsc_snps_to_check (in=in1) remaining_tests (in=in2);
   by snp_id;
   if in1 and in2 then output ucsc_snps_tested;
   else if in1 then output no_test;
   else output no_ucsc;
run;

* Additional 4 SNPs to include;
* 893 Onengut SNPs without a test -- are these in non-expressed genes, or do they have a low MAF?;

/* Check if specific index SNPs are excluded */
data onengut_snp1;
  set onengut_snp_list_tested;
  keep snp_id;
run;

data onengut_snp2;
  set ucsc_snps_tested;
  keep rs_snp_id;
  rename rs_snp_id=snp_id;
run;

data onengut_snp3;
  set no_test;
  keep rs_snp_id;
  rename rs_snp_id=snp_id;
run;


proc sort data=onengut_snp1 nodup;
   by snp_id;
proc sort data=onengut_snp2 nodup;
   by snp_id;
proc sort data=onengut_snp3 nodup;
   by snp_id;
run;

data onengut_snps;
  length list $12.;
  set onengut_snp1 (in=in1) onengut_snp2 (in=in2) onengut_snp3 (in=in3);
  if in1 then list='snp_test';
  if in2 then list='ucsc_test';
  if in3 then list='no_test';
run;

data index_snps;
   set eqtl.supptable1_snp_list;
   keep snp_id index_snp_rs;
run;

proc sort data=onengut_snps nodup;
   by snp_id;
proc sort data=index_snps nodup;
   by snp_id index_snp_rs;
run;

data onengut_snps2index no_index_oops;
   merge onengut_snps (in=in1) index_snps (in=in2);
   by snp_id;
   if in1 and in2 then output onengut_snps2index;
   else if in1 then output no_index_oops;
run;

data index_no_test index_test;
   set onengut_snps2index;
   if list='no_test' then output index_no_test;
   else output index_test;
   keep index_snp_rs;
run;

proc sort data=index_no_test nodup;
   by index_snp_rs;
proc sort data=index_test nodup;
   by index_snp_rs;
run;

data index_snps_testable;
    merge index_test (in=in1) index_no_test;
    by index_snp_rs;
    if in1 then flag_tested=1;
    else flag_tested=0;
run;


/* Stack matching SNPs and make dataset permenant */

data ucsc_snps_tested_2;
  set ucsc_snps_tested;
  rename snp_id=ucsc_snp_id rs_snp_id=snp_id;
run;

data eqtl.onengut_snps_w_tests;
   set onengut_snp_list_tested ucsc_snps_tested_2;
run;







/* The list of non-testable index SNPs are:
SNP_ID		GENE (paper)		PREVIOUS SNP INDEX (paper)
rs10277986	intergenic region	rs4948088 (0.86), rs10231420 (<0.1)
rs34536443	TYK2			rs2304256 (r2<0.1)
rs35667974	IFIH1			rs1990760 (<0.1)	Note: conditional on:rs2111485
rs61839660	IL2RA			rs7090530 (<0.1), rs12251307 (0.61)
rs689		INS			rs7111341 (0.265)

Check if the previous indices are on IC:
rs689 and rs34536443 not on chip!!
rs10277986, rs35667974,rs61839660 have MAF <5% so are filtered out
rs4948088, rs2304256, rs1990760, rs12251307, rs7111341 on IC
rs2304256, rs1990760, rs7111341 on IC after MAF filtering

"SNPs from five regions were excluded on the basis of minor allele frequency, or were missing from the available SNP data."

"Statistically indistinguishable" SNPs for rs689 and rs34536443 also not on IChip
*/

