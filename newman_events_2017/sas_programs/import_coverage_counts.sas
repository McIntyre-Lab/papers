/* Import counts */

libname conesa '!MCLAB/conesa_pacbio/sas_data';
%include '!MCLAB/conesa_pacbio/sas_programs/macros/iterdataset.sas';

data sample_list;
   format sample $4.;
   input sample $;
   datalines; 
   NSC1
   NSC2
   OLD1
   OLD2
   ;
run;


%macro import_counts(sample);

proc import datafile="!MCLAB/conesa_pacbio/alignment_output/coverage_counts_chunks/cvrg_cnts_&sample..csv"
       out=chunk_&sample. dbms=csv replace; guessingrows=100000;
run;

proc import datafile="!MCLAB/conesa_pacbio/alignment_output/coverage_counts_fusions/cvrg_cnts_&sample..csv"
       out=fusion_&sample. dbms=csv replace; guessingrows=100000;
run;

proc import datafile="!MCLAB/conesa_pacbio/alignment_output/coverage_counts_splicing/cvrg_cnts_&sample..csv"
       out=splicing_&sample. dbms=csv replace; guessingrows=686000;
run;

%mend;

%iterdataset(dataset=sample_list, function=%nrstr(%import_counts(&sample);));


data conesa.coverage_fusions;
   set fusion_: ;
run;

data conesa.coverage_exon_chunks;
   set chunk_: ;
run;

data conesa.coverage_splicing;
   set splicing_: ;
run;

