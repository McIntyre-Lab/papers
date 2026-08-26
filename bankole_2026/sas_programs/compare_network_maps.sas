



PROC IMPORT OUT= components_v1
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

PROC IMPORT OUT= components_v2
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/component_map_by_node_atstep2.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data components_v1_4_merge;
set components_v1;
rename component_ID=componentID_v1;
keep jxnhash component_id;
run;

proc freq data=components_v2 noprint;
tables component_id/out=count_compid;
run;


proc sort data=components_v1_4_merge;
by jxnhash;

proc sort data=components_v2;
by jxnhash;

data compare;
merge components_v1_4_merge components_v2;
by jxnhash;
run;

proc freq data=compare noprint;
tables component_id*componentID_v1/out=compare_comp;
run;


proc freq data=compare_comp noprint;
tables component_id/out=count_v2;
tables componentid_v1/out=count_v1;
run;
run;


proc freq data=count_v2;
tables count;
run;

The FREQ Procedure

Frequency Count
COUNT 	Frequency 	Percent 	Cumulative
Frequency 	Cumulative
Percent
1 	61722 	96.15 	61722 	96.15
2 	2013 	3.14 	63735 	99.29
3 	350 	0.55 	64085 	99.83
4 	78 	0.12 	64163 	99.95
5 	27 	0.04 	64190 	100.00
6 	2 	0.00 	64192 	100.00
7 	1 	0.00 	64193 	100.00



proc freq data=count_v1;
tables count;
run;


Frequency Count
COUNT 	Frequency 	Percent 	Cumulative
Frequency 	Cumulative
Percent
1 	49233 	87.34 	49233 	87.34
2 	4995 	8.86 	54228 	96.20
3 	1230 	2.18 	55458 	98.38
4 	524 	0.93 	55982 	99.31
5 	239 	0.42 	56221 	99.74
6 	80 	0.14 	56301 	99.88
7 	22 	0.04 	56323 	99.92
8 	22 	0.04 	56345 	99.96
9 	13 	0.02 	56358 	99.98
10 	7 	0.01 	56365 	99.99
11 	1 	0.00 	56366 	99.99
12 	2 	0.00 	56368 	100.00
14 	2 	0.00 	56370 	100.00


#see how many components are unchanged;

data count_v1a;
set count_v1;
rename count=num_v1_inv2;
drop percent;

data count_v2a;
set count_v2;
rename count=num_v2_inv1;
drop percent;
run;

proc sort data=count_v1a;
by componentid_v1;

proc sort data=compare_comp;
by componentid_v1;
/*56370*/

data frequency_with_counts;
merge compare_comp (in=in1) count_v1a (in=in2);
by componentid_v1;
if in1 and in2;
run;

proc sort data=count_v2a;
by component_id;


proc sort data=frequency_with_counts;
by component_id;
/*64193*/
data frequency_with_counts2;
merge frequency_with_counts (in=in1) count_v2a (in=in2);
by component_id;
if in1 and in2;
run;

proc freq data=frequency_with_counts2;
tables num_v1_inv2*num_v2_inv1;
run;

/*44557 unchanged; 
9352 that are 2 components at step 2 and 1 component in v1 ;
3549 that were 3 components at step 2 and 1 component in v1 ;
2066 that were 4 ";
1166 that were 5 ";
470 that were 6;
148 that were 7 ";
170 that were 8;
112 that were 9;
69 that were 10;
11 that were 11;
24 that were 12;
28 that were 14;
total 61722;





data check23;
set compare;
where componentid_v1=23;
run;
 
 #geneset 23;


