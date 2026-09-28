


PROC IMPORT OUT= mel_gene
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/gene_summary_dmel.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data mel_gene2;
set mel_gene;
dmel_geneid=dmel_fbgn;

PROC IMPORT OUT= mel_dup
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/duplication_vs_splicing/analysis_adp/dmel_Hahn_multigene_classification_with_dup_splicing.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;


data mel_dup2;
set mel_dup;
dmel_geneid=dmel_fbgn;


PROC IMPORT OUT= mel_as
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/summary_files/gene_summary_from_as_analysis_dmel.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data mel_as2;
set mel_as;
rename geneid=dmel_geneid;
run;


PROC IMPORT OUT= mel_anno
            DATAFILE= "/nfshome/mcintyre/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_dmel6_full_annotation.csv"
DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data mel_anno2;
set mel_anno;
rename geneid=dmel_geneid;
run;

/* compare transcript number and gene duplicaiton */

proc sort data=mel_gene2;
by dmel_geneid;

proc sort data=mel_anno2;
by dmel_geneid;

data mel_gene_plus;
merge mel_gene2 mel_anno2;
by dmel_geneid;

num_transcripts= num_ujc -num_ism_ujc ;

run;

proc sort data=mel_dup;
by dmel_fbgn;
run;


data mel_gene_plus_dup nodata;
merge mel_gene2 (in=in1) mel_dup2 (in=in2);
by dmel_geneid;
if in1 then output  mel_gene_plus_dup ;
else if in2 then output nodata;
run;

proc sort data=mel_gene_plus_dup ;
by xscript_class;

proc freq data=mel_gene_plus_dup ;
by xscript_class;
tables flag_gene_de*class_type;
run;

