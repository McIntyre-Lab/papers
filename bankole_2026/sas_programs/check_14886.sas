

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


data mel_14886;
set mel_anno;
where component_id=14886;
jxnhash=dmel6_jxnhash;
run;


data sim_14886;
set sim_anno;
where component_id=14886;
jxnhash=dsim2_jxnhash;
run;


data ser_14886;
set ser_anno;
where component_id=14886;
jxnhash=dser1_jxnhash;
run;



data yak_14886;
set yak_anno;
where component_id=14886;
run;



data san_14886;
set san_anno;
where component_id=14886;
run;



data ser_14886;
set ser_anno;
where component_id=14886;
jxnhash=dser1_jxnhash;
run;



data comp_14886;
set all_componets_w_topology;
where component_id=14886;
run;


PROC IMPORT OUT= ser_start
DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dser1_anno_files/dser11_2_dser1_ujc_xscript_link.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

proc contents data=ser_start;
run;


proc sort data=ser_start;
by jxnhash;

proc sort data= ser_14886;
by jxnhash;

data check_origin;
merge ser_start ser_14886(in=in2);
by jxnhash;
if in2;
run;

proc print data=check_origin;
var transcriptid;
run;



PROC IMPORT OUT= mel_start
DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dmel6_anno_files/dmel650_2_dmel6_ujc_xscript_link.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;



proc sort data=mel_start;
by jxnhash;

proc sort data= mel_14886;
by jxnhash;

data check_mel_origin;
merge mel_start mel_14886(in=in2);
by jxnhash;
if in2;
run;

proc print data=check_mel_origin;
var geneid transcriptid;
run;

/* FBgn0287825 FBtr0087176 */


PROC IMPORT OUT= sim2_start
DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dsim2_anno_files/dsim202_2_dsim2_ujc_xscript_link.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


proc sort data=sim2_start;
by jxnhash;

proc sort data= sim_14886;
by jxnhash;

data check_sim2_origin;
merge sim2_start sim_14886(in=in2);
by jxnhash;
if in2;
run;

proc print data=check_sim2_origin;
var geneid transcriptid;
run;

/* FBgn0196831    FBtr0225455*/

PROC IMPORT OUT= simw_start
DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_supporting_files/fiveSpecies_2_dsim2_anno_files/dsimWXD_2_dsim2_ujc_xscript_link.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


proc sort data=simw_start;
by jxnhash;

proc sort data= sim_14886;
by jxnhash;

data check_simw_origin;
merge simw_start sim_14886(in=in2);
by jxnhash;
if in2;
run;

proc print data=check_simw_origin;
var geneid transcriptid;
run;

/*  FBgn0196831    maker-2R-exonerate_est2genome-gene-48.29-mRNA-2*/
