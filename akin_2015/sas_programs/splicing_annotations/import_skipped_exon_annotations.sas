/******** IMPORT SKIPPED EXON ANNOTATIONS  **********/


libname splice '!MCLAB/junction_annotations/sas_data';

/* Import skipped exon annotations */

    data WORK.EXON_SKIPPING_ANNOT    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exon_skipping_annot.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat junction_id $79. ;
       informat flag_exonskip best32. ;
       informat num_skipped_exons best32. ;
       informat cat_skipped_exons $7489. ;
       format junction_id $79. ;
       format flag_exonskip best12. ;
       format num_skipped_exons best12. ;
       format cat_skipped_exons $7489. ;
    input
                junction_id $
                flag_exonskip
                num_skipped_exons
                cat_skipped_exons $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


/*     data WORK.SKIPPED_EXON_LIST    ; */
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
/*     infile
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_skipped_exons_list.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat junction_id $79. ;
        informat skipped_exon_id $39. ;
        informat flag_exonskip best32. ;
        format junction_id $79. ;
        format skipped_exon_id $39. ;
        format flag_exonskip best12. ;
     input
                 junction_id $
                 skipped_exon_id $
                 flag_exonskip
     ;
     if _ERROR_ then call symputx('_EFIERR_',1); */ /* set ERROR detection macro variable */
/*     run;  */


/* Make permenant */

data splice.EXON_SKIPPING_ANNOT;
  set EXON_SKIPPING_ANNOT;
  drop cat_skipped_exons; *don't need, but can reimport if necessary;
run;

/* if I need this later, I can reimport it
data splice.skipped_exon_list;
  set skipped_exon_list;
run; */

proc datasets noprint;
  delete EXON_SKIPPING_ANNOT skipped_exon_list;
run;
quit;
