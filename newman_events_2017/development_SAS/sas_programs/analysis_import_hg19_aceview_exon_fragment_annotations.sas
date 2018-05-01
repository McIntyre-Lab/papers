ods listing; ods html close;
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Import exon info */

    data WORK.EXONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!PATCON/useful_human_data/aceview_hg19/created_files/hg19_aceview_exons.csv' delimiter = ','
MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chrom $5. ;
       informat start best32. ;
       informat stop best32. ;
       informat strand $1. ;
       informat exon_id $39. ;
       informat transcript_id $1522. ;
       informat gene_id $36. ;
       format chrom $5. ;
       format start best12. ;
       format stop best12. ;
       format strand $1. ;
       format exon_id $39. ;
       format transcript_id $1522. ;
       format gene_id $36. ;
    input
                chrom $
                start
                stop
                strand $
                exon_id $
                transcript_id $
                gene_id $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


/* Import exon fragment info */

    data WORK.EXON_FRAGMENTS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!PATCON/useful_human_data/aceview_hg19/created_files/hg19_aceview_exon_fragments.csv' delimiter =
 ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chrom $5. ;
       informat group_start best32. ;
       informat group_stop best32. ;
       informat fragment_start best32. ;
       informat fragment_end best32. ;
       informat exon_id $805. ;
       format chrom $5. ;
       format group_start best12. ;
       format group_stop best12. ;
       format fragment_start best12. ;
       format fragment_end best12. ;
       format exon_id $805. ;
    input
                chrom $
                group_start
                group_stop
                fragment_start
                fragment_end
                exon_id $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;



data hg19.hg19_aceview_exons;
   set exons;
run;

data hg19.hg19_aceview_exon_fragments;
   set exon_fragments;
run;


