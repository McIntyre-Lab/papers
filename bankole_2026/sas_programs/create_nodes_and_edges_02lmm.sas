

%macro import_genome (genome);

/* import full anno files */
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

proc contents node_list;
run;


/*length of geneid varies which one is longest think it is ser*/

%macro import_links2 (anno, genome);

PROC IMPORT OUT= WORK.&anno._2_&genome
            DATAFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_&anno._anno_files/fivespecies_&anno._2_&genome._noGeneID_ujc_xscript_link.csv" 
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


%import_links2 (dmel6,dsan1);

%import_links2 (dmel6,dser1);
%import_links2 (dmel6,dyak2);
%import_links2 (dmel6,dsim2);

%import_links2 (dsim2,dsan1);
%import_links2 (dsim2,dser1);
%import_links2 (dsim2,dyak2);
%import_links2 (dsim2,dmel6);

%import_links2 (dsan1,dmel6);
%import_links2 (dsan1,dser1);
%import_links2 (dsan1,dyak2);
%import_links2 (dsan1,dsim2);

%import_links2 (dyak2,dsan1);
%import_links2 (dyak2,dser1);
%import_links2 (dyak2,dmel6);
%import_links2 (dyak2,dsim2);

%import_links2 (dser1,dsan1);
%import_links2 (dser1,dmel6);
%import_links2 (dser1,dyak2);
%import_links2 (dser1,dsim2);

data links_big;
set 

dmel6_2_dsan1_id

dmel6_2_dser1_id
 dmel6_2_dyak2_id
 dmel6_2_dsim2_id

 dsim2_2_dsan1_id
 dsim2_2_dser1_id
 dsim2_2_dyak2_id
 dsim2_2_dmel6_id

 dsan1_2_dmel6_id
 dsan1_2_dser1_id
 dsan1_2_dyak2_id
 dsan1_2_dsim2_id

 dyak2_2_dsan1_id
 dyak2_2_dser1_id
 dyak2_2_dmel6_id
 dyak2_2_dsim2_id

 dser1_2_dsan1_id
 dser1_2_dmel6_id
 dser1_2_dyak2_id
dser1_2_dsim2_id

;
run;

PROC EXPORT DATA= links
            OUTFILE= "!MCLAB/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files/edges.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
