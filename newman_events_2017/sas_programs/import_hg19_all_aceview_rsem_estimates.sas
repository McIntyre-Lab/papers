ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname con '!PATCON/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
filename mymacros '!MCLAB/event_analysis/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Import RSEM results for the full hg19 Aceview transcriptome */

* Get list of libraries, as I am using this for sample ID here;

%include "!MCLAB/event_analysis/sas_programs/mclib_SAS/iterdataset.sas";

data libraries;
  set con.design_by_subject_new;
  keep library;
run;


proc datasets noprint;
  delete rsem_: SL:;
run; quit;


%macro importRSEM(library);

    data WORK.&library._est   ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
"/mnt/store/event_sandbox/rsem_output/aceview_hg19_all/_&library..isoforms.results"
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat transcript_id $113. ;
       informat gene_id $36. ;
       informat length best32. ;
       informat effective_length best32. ;
       informat expected_count best32. ;
       informat TPM best32. ;
       informat FPKM best32. ;
       informat IsoPct best32. ;
       format transcript_id $113. ;
       format gene_id $36. ;
       format length best12. ;
       format effective_length best12. ;
       format expected_count best12. ;
       format TPM best12. ;
       format FPKM best12. ;
       format IsoPct best12. ;
    input
                transcript_id $
                gene_id $
                length
                effective_length
                expected_count
                TPM
                FPKM
                IsoPct
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


*trim data -- only want TPM estimate and transcript_id;

data rsem_&library.;
  length library $8.;
  set &library._est;
  library="&library.";
  keep transcript_id library TPM;
run;

%mend;
   %iterdataset(dataset=libraries, function=%nrstr(%importRSEM(&library);));

data all_xscript_est;
   set rsem_SL: ;
run;

/* Update IDs */

data xs2fasta;
   set event.hg19_aceview_xs2gene_fasta_index;
   keep transcript_id transcript_fasta_id;
   rename transcript_id=aceview_xs_id transcript_fasta_id=transcript_id;
run;

proc sort data=all_xscript_est;
   by transcript_id;
proc sort data=xs2fasta;
   by transcript_id;
run;

data all_xscript_est2;
  merge all_xscript_est (in=in1) xs2fasta (in=in2);
   by transcript_id;
  if in1 and in2;
run;

/* Make permenant */

data eventloc.hg19_rsem_all_xscripts;
   set all_xscript_est2;
   drop transcript_id;
   rename aceview_xs_id=transcript_id;
run;


SL30351
