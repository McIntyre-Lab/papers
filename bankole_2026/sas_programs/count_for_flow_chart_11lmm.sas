

PROC IMPORT OUT= mel_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dmel6_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


PROC IMPORT OUT= sim_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dsim2_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


PROC IMPORT OUT= yak_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dyak2_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

PROC IMPORT OUT= san_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dsan1_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

PROC IMPORT OUT= ser_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dser1_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;



PROC IMPORT OUT= components
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_network_files/component_map_by_node.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
run;


PROC IMPORT OUT= component_summary
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/component_summary.csv"
	DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


/* the origin of the jxnhash relative to each species!*/

data mel_origin;
set mel_anno;
where flag_dmel650_2_dmel6_ujc=1;
rename dmel6_jxnhash=jxnhash;
keep dmel6_jxnhash;
run;

proc sort data=mel_origin ;
by jxnhash;
run;

proc sort data=components;
by jxnhash;

data component_mel_origin;
merge mel_origin (in=in1) components (in=in6);
by jxnhash;
if in1;
run;

proc sort data= component_mel_origin /*nodupkey*/;
by component_id;
run;

data sim_origin;
set sim_anno;
where flag_dsim202_2_dsim2_ujc=1 or flag_dsimWXD_2_dsim2_ujc=1;
rename dsim2_jxnhash=jxnhash;
keep dsim2_jxnhash;
run;

proc sort data=sim_origin ;
by jxnhash;
run;


data component_sim_origin;
merge sim_origin (in=in1) components (in=in6);
by jxnhash;
if in1;

proc sort data= component_sim_origin nodupkey;
by component_id;
run;

data yak_origin;
set yak_anno;
where flag_dyak21_2_dyak2_ujc=1;
rename dyak2_jxnhash=jxnhash;
keep dyak2_jxnhash;
run;

proc sort data=yak_origin ;
by jxnhash;
run;

data component_yak_origin;
merge yak_origin (in=in1) components (in=in6);
by jxnhash;
if in1;

proc sort data= component_yak_origin nodupkey;
by component_id;
run;


data san_origin;
set san_anno;
where flag_dsan11_2_dsan1_ujc=1;
rename dsan1_jxnhash=jxnhash;
keep dsan1_jxnhash;
run;

proc sort data=san_origin ;
by jxnhash;
run;

data component_san_origin;
merge san_origin (in=in1) components (in=in6);
by jxnhash;
if in1;

proc sort data= component_san_origin nodupkey;
by component_id;
run;

data ser_origin;
set ser_anno;
where flag_dser11_2_dser1_ujc=1;
rename dser1_jxnhash=jxnhash;
keep dser1_jxnhash;
run;

proc sort data=ser_origin ;
by jxnhash;
run;

data component_ser_origin;
merge ser_origin (in=in1) components (in=in6);
by jxnhash;
if in1;

proc sort data= component_ser_origin nodupkey;
by component_id;
run;


data component_origin ;
merge component_mel_origin (in=in1) component_sim_origin(in=in2) component_yak_origin(in=in3) component_san_origin(in=in4) component_ser_origin(in=in5);
by component_id;

num_species_origin=in1+in2+in3+in4+in5;

origin_pattern=(cats(in1,in2,in3,in4,in5));

if in1 and in2 and in3 and in4 and in5 then origin="1_all_5_species";
else if in1 and in2 and in3 and in4 then origin="2_the4";
else if in1 and in2 and num_species_origin=2 then origin="3_mel_sim";
else if in3 and in4 and num_species_origin=2 then origin="3_yak_san";
else if in1 and num_species_origin=1 then origin="4_melonly";
else if in2 and num_species_origin=1 then origin="4_simonly";
else if in3 and num_species_origin=1 then origin="4_yakonly";
else if in4 and num_species_origin=1 then origin="4_sanonly";
else if in5 and num_species_origin=1 then origin="4_seronly";
componentid=component_id;
keep componentid origin origin_pattern;
run;

proc sort data=component_origin nodupkey;
by componentid origin origin_pattern;
run;

/*figure out which species are associated with each of the components in the annotation*/

proc freq data=components noprint;
tables component_id*source/out=count_comps_by_source;
run;


/*56,370 components*/


proc freq data=count_comps_by_source noprint;
tables component_id/out=count_source;
run;

proc freq data=count_comps_by_source ;
tables source/out=count_comps2;
run;

proc freq data=count_source;
tables count;
run;

proc transpose data=count_comps_by_source out=flip;
by component_id;
var count;
id source;
run;

/*components linked to dmel*/

proc freq data=flip;
where dmel6>0;
tables dsim2 dsan1 dyak2 dser1;
run;

%macro components_by_species (species);

data components_4_&species._merge;
set components;
&species._jxnhash=jxnhash;
where source="&species";
drop geneID source jxnhash;
run;

proc freq data=components_4_&species._merge noprint;
tables component_id/out=&species._num_jxnhash;
run;

data components_4_&species._merge1 ;
set &species._num_jxnhash;
rename count=num_&species._jxnhash;
drop percent;
run;
%mend;

%components_by_species (dmel6); /*75986 jxnhash 49124 components */
%components_by_species (dsim2); /*77281 jxnhash 49240*/
%components_by_species (dser1); /*72272 jxnhash 44898*/
%components_by_species (dyak2); /*77040 jxnhash 49088*/
%components_by_species (dsan1); /*77023 jxnhash 49045*/

data componets_by_species;
length comp_anno_in $20;
merge components_4_dmel6_merge1 (in=in1)
components_4_dsim2_merge1 (in=in2)
components_4_dyak2_merge1 (in=in3)
components_4_dsan1_merge1 (in=in4)
components_4_dser1_merge1 (in=in5)
;
by component_id;
num_species=in1+in2+in3+in4+in5;
if in1 and in2 and in3 and in4 and in5 then comp_anno_in="1_all_5_species";
else if in1 and in2 and in3 and in4 then comp_anno_in="2_the4";
else if in1 and in2 and num_species=2 then comp_anno_in="3_mel_sim";
else if in3 and in4 and num_species=2 then comp_anno_in="3_yak_san";
else if in1 and num_species=1 then comp_anno_in="4_melonly";
else if in2 and num_species=1 then comp_anno_in="4_simonly";
else if in3 and num_species=1 then comp_anno_in="4_yakonly";
else if in4 and num_species=1 then comp_anno_in="4_sanonly";
else if in5 and num_species=1 then comp_anno_in="4_seronly";
else comp_anno_in="6_other";
componentid=component_id;
run;

/*56370*/

proc freq data=componets_by_species;
tables comp_anno_in;
run;

proc sort data=component_summary;
by componentid;


data all_components2;
length expressed_in $20;
merge component_origin componets_by_species (in=in1) component_summary (in=in2);
by componentid;
if flag_expressed=0 then expressed_in="7_not_expressed";
else if sumReads_T_dmel >0 and sumReads_T_dsim>0 and sumReads_T_dyak >0 
	and sumReads_T_dsan >0 and sumReads_T_dser >0 then expressed_in="1_all_5";

else if sumReads_T_dmel >0 and sumReads_T_dsim>0 and sumReads_T_dyak >0 
	and sumReads_T_dsan >0 then expressed_in="2_the4";

else if sumReads_T_dmel >0 and sumReads_T_dsim>0 and sumReads_T_dyak =0 
	and sumReads_T_dsan =0 and sumReads_T_dser =0 then expressed_in="3_mel_sim";

else if sumReads_T_dmel =0 and sumReads_T_dsim=0 and sumReads_T_dyak >0 
	and sumReads_T_dsan >0 and sumReads_T_dser =0  then expressed_in="3_yak_san";

else if sumReads_T_dmel >0 and sumReads_T_dsim=0 and sumReads_T_dyak =0 
	and sumReads_T_dsan =0 and sumReads_T_dser =0 then expressed_in="4_mel_only";

else if sumReads_T_dmel =0 and sumReads_T_dsim>0 and sumReads_T_dyak =0 
	and sumReads_T_dsan =0 and sumReads_T_dser =0  then expressed_in="4_sim_only";

else if sumReads_T_dmel =0 and sumReads_T_dsim=0 and sumReads_T_dyak =0 
	and sumReads_T_dsan >0 and sumReads_T_dser =0  then expressed_in="4_san_only";

else if sumReads_T_dmel =0 and sumReads_T_dsim=0 and sumReads_T_dyak >0 
	and sumReads_T_dsan =0 and sumReads_T_dser =0  then expressed_in="4_yak_only";

else if sumReads_T_dmel =0 and sumReads_T_dsim=0 and sumReads_T_dyak =0 
	and sumReads_T_dsan =0 and sumReads_T_dser >0  then expressed_in="4_ser_only";
else expressed_in="6_other";

if origin ="" then origin="6_other";
if origin_pattern="" then origin_pattern="00000";

anno_in=substr(comp_anno_in,1,2);
expr_in=substr(expressed_in,1,2);
run;


proc freq data=all_components2;
where comp_anno_in="1_all_5";
*tables expressed_in;
tables expr_in;
run;

data component_info;
set all_components2;
keep component_id origin origin_pattern comp_anno_in expressed_in flag_expressed ;
run;

proc freq data=all_components2 noprint;
tables origin*comp_anno_in*flag_expressed/out=comp_orgn_anno_xpr;
run;
quit;

proc export data = comp_orgn_anno_xpr
outfile = "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/count_origin_anno_xpress.csv"
dbms = csv replace ;
run;

proc freq data=all_components2 noprint;
tables origin*comp_anno_in*expressed_in/out=comp_orgn_anno_xpr_byspecies;
run;

proc export data = comp_orgn_anno_xpr_byspecies
outfile = "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/count_origin_anno_xpress_byspecies.csv"
dbms = csv replace ;
run;


/* after this is checks and various tabulations to look for interesting transcripts not used in results component_info should be replaced by all_components_w_topology*/


proc export data = component_info
outfile = "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/component_info.csv"
dbms = csv replace ;
run;






