libname splicing '/mnt/data/splicing/';
libname splice '/mnt/data/splice/';

/* Adding  in event annotations */

data splicing_annotations;
  set splice.splicing_events_annotations;
  keep event_id gene_id flag_junction_annotated flag_intron_retention flag_exonskip flag_alt_donor flag_alt_acceptor;
run;


/* sort data */
proc sort data=splicing.results_by_splicing_w_fdr;
    by event_id;
run;


proc sort data=splicing_annotations;
    by event_id;
run;


/* Merge */

data results_w_annot no_annot;
    merge splicing.results_by_splicing_w_fdr (in=in1) splicing_annotations (in=in2);
    by event_id;
    if in1 and in2 then output results_w_annot;
    else if in1 then output no_annot; *should be 0 obs;
    else do;  *anything with 0 coverage in all three cell types should go here;
     mean_depth_CD19=0;
     mean_depth_CD4=0;
     mean_depth_CD8=0;
     count_CD19=0;
     flag_CD19_on=0;
     count_CD4=0;
     flag_CD4_on=0;
     count_CD8=0;
     flag_CD8_on=0;
     flag_all_on=0;
     output results_w_annot;
     end;
run;

/* Make permenant */

data splicing.splicing_results_w_annot;
   set results_w_annot;
run;



