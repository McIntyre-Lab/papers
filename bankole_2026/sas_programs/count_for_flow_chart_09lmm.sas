

/*all_components is from create_component_annotation.sas*/

/*figure out which species are associated with each of the components in the annotation*/

proc freq data=all_components noprint;
tables component_id*origin/out=count_comps_by_source;
run;


/*56,370 components*/



proc freq data=count_comps_by_source ;
tables origin/out=count_comps2;
run;

proc freq data=count_source;
tables count;
run;
/*
dmel6 49124 20.35 49124 20.35 
dsan1 49045 20.32 98169 40.67 
dser1 44898 18.60 143067 59.27 
dsim2 49240 20.40 192307 79.66 
dyak2 49088 20.34 241395 100.00 */




/*components linked to dmel*/


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
if in1 and in2 and in3 and in4 and in5 then comp_anno_in="1_all_5";
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
run;

proc contents data=component_summary;
run;

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
run;


proc freq data=all_components2;
tables flag_expressed;
run;

data component_info;
set all_components;
keep component_id origin origin_pattern comp_anno_in expressed_in ;
run;

proc freq data=all_components noprint;
tables origin*comp_anno_in*flag_expressed/out=comp_orgn_anno_xpr;
run;
quit;

proc export data = comp_orgn_anno_xpr
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\comp_orgn_anno_xpr.csv"
dbms = csv replace ;
run;

proc freq data=all_components noprint;
tables origin*comp_anno_in*expressed_in/out=comp_orgn_anno_xpr_byspecies;
run;

proc export data = comp_orgn_anno_xpr_byspecies
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\comp_orgn_anno_xpr_byspecies.csv"
dbms = csv replace ;
run;


/* after this is checks and various tabulations to look for interesting transcripts not used in results  component_info should be replaced by all_components_w_topology*/

/*
proc export data = component_info
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\submission\supplementary\fiveSpecies_annotations\link_files\component_info.csv"
dbms = csv replace ;
run;*/


proc freq data=all_components noprint;
where flag_expressed=1;
tables comp_anno_in*expressed_in/out=comp_anno_xpr;
run;


proc export data = comp_anno_xpr
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\comp_anno_xpr.csv"
dbms = csv replace ;
run;


proc freq data=all_components;
where comp_anno_in="1_all_5";
tables expressed_in;
run;

proc freq data=all_components;
where expressed_in="1_all_5";
tables comp_anno_in;
run;



proc export data = cnt_mel_sim_sex
outfile = "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/cnt_mel_sim_sex.csv"
dbms = csv replace ;
run;



proc freq data=all_components noprint;
where dmel_num_genes >0 and dsim_num_genes >0 and dser_num_genes>0;
tables origin*comp_anno_in*dmel_sexclass*dsim_sexclass*dser_sexclass/out=cnt_mel_sim_ser_sex;
run;



proc export data = cnt_mel_sim_ser_sex
outfile = "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/cnt_mel_sim_ser_sex.csv"
dbms = csv replace ;
run;


data check;
set all_components;

where flag_expressed=1 and origin="4_zseronly" and dmel_sexclass=dsim_sexclass and dmel_sexclass ne "DE_neither";
run;

proc print data=check;
var component_id;
run;


proc freq data=all_components;
where comp_anno_in="1_all
tables flag_expressed;
run;

proc freq data=all_components;
tables dmel_sexclass*dsim_sexclass*dser_sexclass/out=count_de_bycomp;
run;

data check_3species_sexbias;
set all_components;
where dmel_sexclass=dsim_sexclass and dmel_sexclass=dser_sexclass;
if dmel_sexclass="DE_neither" then delete;
if dmel_sexclass="" then delete;
run;

proc sort data=check_3species_sexbias;
by dmel_geneid_concat;
run;

data pull_id;
set all_components;
fbgn=substr(dmel_geneid_concat, 1,11);
run;

proc sort data=mel_go;
by fbgn;
run;

proc sort data=pull_id;
by fbgn;

data peek;
merge  pull_id (in=in1) mel_go;
by fbgn;
if in1;
run;

data peek1;
set peek;
where dmel_sexclass=dsim_sexclass and dmel_sexclass=dser_sexclass;
if dmel_sexclass="DE_neither" then delete;
if dmel_sexclass="" then delete;
run;

/*hoip in this list*/

/*no oops*/


proc freq data=all_components noprint;
where flag_expressed=1;
tables comp_anno_in*expressed_in/out=count_combos;
run;





proc freq data=all_components ;
where flag_expressed=1;
tables expressed_in;
run;

proc freq data=all_components ;
where flag_expressed=1;
tables dmel_sexclass*dsim_sexclass;
run;

data check;
set all_components;
 where  dmel_sexclass = dsim_sexclass and  dmel_sexclass="DE_both";
run;


proc freq data=all_components ;
where comp_anno_in="4_seronly";
tables expressed_in;
run;

proc freq data=all_components ;
where expressed_in="5_ser_only";
tables comp_anno_in;
run;


proc export data = count_combos 
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\pattern_components_anno_expressed.csv"
dbms = csv replace ;
run;

proc freq data=all_components noprint;
tables comp_anno_in*flag_expressed/out=count_expressed;
run;

data count_expressed1;
set count_expressed;
if flag_expressed=0 then expressed="not_expressed";
if flag_expressed=1 then expressed="expressed";

proc export data = count_expressed1 
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\counts_components_anno_expressed.csv"
dbms = csv replace ;
run;


proc export data = count_de_bycomp 
outfile = "Z:\SHARE\McIntyre_Lab\sex_specific_splicing\Tables\count_de_bycomp.csv"
dbms = csv replace ;
run;


data count_nomel;
set count_mel_collapsed;
length expressed_in $20;
count_expressed=_freq_-sum_dmel6_jxnhash;
if after="all_5_species" then expressed_in="not_mel";
else if after="mel_sim_yak_san" then expressed_in="sim_yak_san";
else if after="mel_sim" then expressed_in="sim";
else expressed_in="other_nomel";
keep before  expressed_in count_expressed;
run;

data count_mel_collapsed1;
set count_mel_collapsed;
 count_expressed=sum_dmel6_jxnhash;
expressed_in=after;
keep before expressed_in count_expressed;
run;

data more_links;
set count_nomel count_mel_collapsed1;
run;


proc export data = count_mel_collapsed 
outfile = "c:\a1stuff\counts_on_dmel6.csv"
dbms = csv replace ;
run;


proc export data = more_links
outfile = "c:\a1stuff\more_links.csv"
dbms = csv replace ;
run;
