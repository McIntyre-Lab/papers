
ods listing; ods html close;
libname con '!PATCON/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* Import MISO results for T1D data and summarize to gene-level */

/* Import cell-type MISO summaries */


%macro sumMISO(cell);

proc import datafile="/mnt/store/miso_sandbox/misopy-0.5.3/miso_summary_SE_&cell./summary/miso_output_SE_&cell..miso_summary"
    out=miso_summary_&cell. dbms=tab replace; guessingrows=40000;
run;
%mend;

%sumMISO(CD4);
%sumMISO(CD8);
%sumMISO(CD19);


/* Import MISO comparisons */

%macro compMISO(cell1,cell2);

proc import datafile="/mnt/store/miso_sandbox/misopy-0.5.3/miso_comparisons/miso_output_SE_&cell1._vs_miso_output_SE_&cell2./bayes-factors/miso_output_SE_&cell1._vs_miso_output_SE_&cell2..miso_bf"
    out=miso_&cell1._&cell2. dbms=tab replace; guessingrows=40000;
run;
%mend;

%compMISO(CD4,CD8);
%compMISO(CD4,CD19);
%compMISO(CD8,CD19);

/* Import event 2 gene index for MISO */

proc import datafile="/mnt/store/miso_sandbox/hg19/hg19/SE.hg19.gff3_to_ensGene2.txt"
   out=miso_event2gene dbms=tab replace; guessingrows=40000;
run;



/* Make permenant */

data event.hg19_miso_event2gene;
   set miso_event2gene;
run;


data event.t1d_miso_summary_cd4;
   set miso_summary_cd4;
run;

data event.t1d_miso_summary_cd8;
   set miso_summary_cd8;
run;

data event.t1d_miso_summary_cd19;
   set miso_summary_cd19;
run;

data event.t1d_miso_compare_cd4cd8;
   set miso_cd4_cd8;
run;

data event.t1d_miso_compare_cd4cd19;
   set miso_cd4_cd19;
run;

data event.t1d_miso_compare_cd8cd19;
   set miso_cd8_cd19;
run;



