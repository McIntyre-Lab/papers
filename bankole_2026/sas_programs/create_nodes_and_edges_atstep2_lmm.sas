

%macro import_genome (genome);

PROC IMPORT OUT= WORK.&genome
            DATAFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fivespecies_&genome._full_annotation_w_dataFlag_gffCompMM.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;


data &genome._id;
set &genome;
source="&genome";
rename &genome._jxnhash=jxnhash;
keep source &genome._jxnhash geneid ;
run;

%mend;

%import_genome(dmel6);
%import_genome(dsim2);
%import_genome(dser1);
%import_genome(dyak2);
%import_genome(dsan1);

data node_list;
set dser1_id dmel6_id dsim2_id dsan1_id dyak2_id ;
run;

PROC EXPORT DATA= node_list
            OUTFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files/node_list.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/*length of geneid varies which one is longest think it is ser*/

%macro import_links2 (genome,anno);

PROC IMPORT OUT= WORK.&anno._2_&genome
            DATAFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_&genome._anno_files/&anno._2_&genome._noGeneID_ujc_xscript_link.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;


data &anno._2_&genome._id;
set &anno._2_&genome;
rename transcriptid=source_jxnHash
jxnHash=target_jxnHash;
source="&anno";
target="&genome";
drop geneid; 
run;

%mend;


%import_links2 (dmel6,dsan11_2_dsan1_ujc);

%import_links2 (dmel6,dser11_2_dser1_ujc);
%import_links2 (dmel6,dyak21_2_dyak2_ujc);
%import_links2 (dmel6,dsim202_2_dsim2_ujc);
%import_links2 (dmel6,dsimWXD_2_dsim2_ujc);


%import_links2 (dsim2,dsan11_2_dsan1_ujc);
%import_links2 (dsim2,dser11_2_dser1_ujc);
%import_links2 (dsim2,dyak21_2_dyak2_ujc);
%import_links2 (dsim2,dmel650_2_dmel6_ujc);

%import_links2 (dsan1,dmel650_2_dmel6_ujc);
%import_links2 (dsan1,dser11_2_dser1_ujc);
%import_links2 (dsan1,dyak21_2_dyak2_ujc);
%import_links2 (dsan1,dsim202_2_dsim2_ujc);
%import_links2 (dsan1,dsimWXD_2_dsim2_ujc);

%import_links2 (dyak2,dsan11_2_dsan1_ujc);
%import_links2 (dyak2,dser11_2_dser1_ujc);
%import_links2 (dyak2,dmel650_2_dmel6_ujc);
%import_links2 (dyak2,dsim202_2_dsim2_ujc);
%import_links2 (dyak2,dsimWXD_2_dsim2_ujc);

%import_links2 (dser1,dsan11_2_dsan1_ujc);
%import_links2 (dser1,dmel650_2_dmel6_ujc);
%import_links2 (dser1,dyak21_2_dyak2_ujc);
%import_links2 (dser1,dsim202_2_dsim2_ujc);
%import_links2 (dser1,dsimWXD_2_dsim2_ujc);

data links;
set 




dmel650_2_dmel6_ujc_2_dsan1_id
dmel650_2_dmel6_ujc_2_dser1_id
dmel650_2_dmel6_ujc_2_dyak2_id
dmel650_2_dmel6_ujc_2_dsim2_id


 dsim202_2_dsim2_ujc_2_dsan1_id
 dsim202_2_dsim2_ujc_2_dser1_id
 dsim202_2_dsim2_ujc_2_dyak2_id
 dsim202_2_dsim2_ujc_2_dmel6_id

 dsimWXD_2_dsim2_ujc_2_dsan1_id
 dsimWXD_2_dsim2_ujc_2_dser1_id
 dsimWXD_2_dsim2_ujc_2_dyak2_id
 dsimWXD_2_dsim2_ujc_2_dmel6_id

 dsan11_2_dsan1_ujc_2_dmel6_id
dsan11_2_dsan1_ujc_2_dser1_id
dsan11_2_dsan1_ujc_2_dyak2_id
dsan11_2_dsan1_ujc_2_dsim2_id


 dyak21_2_dyak2_ujc_2_dsan1_id
 dyak21_2_dyak2_ujc_2_dser1_id
 dyak21_2_dyak2_ujc_2_dmel6_id
 dyak21_2_dyak2_ujc_2_dsim2_id


 dser11_2_dser1_ujc_2_dsan1_id
 dser11_2_dser1_ujc_2_dmel6_id
 dser11_2_dser1_ujc_2_dyak2_id
 dser11_2_dser1_ujc_2_dsim2_id


;
run;



PROC EXPORT DATA= links
            OUTFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files/edges_step2.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc contents data=links;
run;

#compare to links from create_nodes_and_edges_02lmm.sas

proc sort data =links;
by source_jxnhash target_jxnhash;

proc sort data=links_big;
by source_jxnhash target_jxnhash;


data instep2 overlap added_by_remap;
merge links_big(in=in1) links(in=in2);
by source_jxnhash target_jxnhash;
if in1 and in2 then output overlap;
else if in1 then output added_by_remap;
else if in2 then output instep2;
run;


data with_extra_edges;
set links_big instep2;
run;


