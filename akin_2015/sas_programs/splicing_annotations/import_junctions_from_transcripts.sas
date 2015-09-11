/******** IMPORT TRANSCRIPT-ANNOTATED JUNCTIONS  **********/

libname splice '!MCLAB/junction_annotations/sas_data';

/* Import transcript-annotated junctions */

    data WORK.TRANSCRIPT_JUNCTIONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_transcript_junctions.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat junction_id $79. ;
       informat junc_coords $24. ;
       informat transcript_id $53. ;
       informat gene_id $36. ;
       format junction_id $79. ;
       format junc_coords $24. ;
       format transcript_id $53. ;
       format gene_id $36. ;
    input
                junction_id $
                junc_coords $
                transcript_id $
                gene_id $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Format transcript junctions */

data xscript_junc_formatted;
   set transcript_junctions;
   drop junc_coords; * Had this in here for multigene junctions, but going to do this in a later program;
run;


/* Cat together transcripts by junction */

proc sort data=xscript_junc_formatted;
    by junction_id transcript_id;
run;

/* get counts first */

proc freq noprint data=xscript_junc_formatted;
   tables junction_id / out=junc_count;
run;

proc sort data=junc_count;
  by descending count;
run;


*max=55 transcripts per junction;


data junctions_cat_xscript; 
  array xscripts[55] $ 53;

  retain xscripts1-xscripts55;

  set xscript_junc_formatted;
  by junction_id;
  
  if first.junction_id then do;
     call missing(of xscripts1-xscripts55);
     records = 0;
  end;

  records + 1;
  xscripts[records]=transcript_id;
  if last.junction_id then output;
run;

  *clean up the output file;

data junctions_cat_xscript2;
  set junctions_cat_xscript;
  length cat_xscript $ 2970;
  rename records= num_transcripts;
         cat_xscript= catx("|", OF xscripts1-xscripts55);
  drop xscripts1-xscripts55 transcript_id;
  rename cat_xscript=transcript_id;
  run;


/* Make permenant */

data splice.xscript_junctions;
   set junctions_cat_xscript2;
run;

/* remove temp datasets */

proc datasets noprint;
  delete TRANSCRIPT_JUNCTIONS xscript_junc_formatted junc_count
         junctions_cat_xscript junctions_cat_xscript2;
run;
quit;

