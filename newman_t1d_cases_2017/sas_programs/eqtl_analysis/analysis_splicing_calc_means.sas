/*** IMPORT JUNCTION DATA ***/

/* import libraries */

libname mysas '/ufrc/concannon/share/jnewman/sas_analysis/sas_data2';

/**********************************************************************************************************/

/* DEPTH */
data jnc_data_for_tpose_&sysparm.;
   set mysas.counts_by_splicing_w_flags_&sysparm.;
   keep name event_id depth;
run;

proc sort data=jnc_data_for_tpose_&sysparm.;
    by event_id name;
run;

proc transpose data=jnc_data_for_tpose_&sysparm. out=flipped_depth_&sysparm.;
   by event_id;
   var depth;
   id name;
run;

data flipped_depth_2_&sysparm.;
   set flipped_depth_&sysparm.;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   drop _NAME_;
   run ;

/* untranspose! */

proc transpose name=Name data=flipped_depth_2_&sysparm. out=splice_depth_untpose_&sysparm.(rename=(col1=depth));
     by event_id;
     var _2009_PC_:;
run;

data splice_depth_untpose_2_&sysparm.;
   set splice_depth_untpose_&sysparm.;
   Name=tranwrd(Name, "_2009_PC_", "2009-PC-");
run;



data splicing_counts_w_zeros_&sysparm.;
set splice_depth_untpose_2_&sysparm.;
run;


/* make perm */

data mysas.splicing_counts_w_zeros_&sysparm.;
  set splicing_counts_w_zeros_&sysparm.;
run;

/* Add in cell type */

data cell_types_&sysparm.;
  set mysas.counts_by_splicing_w_flags_&sysparm.;
  keep Name cell_type;
run;

proc sort data=cell_types_&sysparm. nodup;
   by Name;
run;

proc sort data=splicing_counts_w_zeros_&sysparm.;
   by Name;
run;

data splicing_counts_w_cell_&sysparm. oops1 oops2;
   merge cell_types_&sysparm. (in=in1) splicing_counts_w_zeros_&sysparm. (in=in2);
   by Name;
   if in1 and in2 then output splicing_counts_w_cell_&sysparm.;
   else if in1 then output oops1;
   else output oops2;
run;


/* Get group means (APN, log_apn) */

proc sort data=splicing_counts_w_cell_&sysparm.;
   by event_id cell_type;
run;

proc means data=splicing_counts_w_cell_&sysparm. noprint;
   by event_id cell_type;
   var depth;
   output out=jnc_counts_means_&sysparm. mean=;
run;

*split on treatment, rename variables and merge;

data jnc_CD19_mean_&sysparm. jnc_CD4_mean_&sysparm. jnc_CD8_mean_&sysparm.;
    set jnc_counts_means_&sysparm.;
    if cell_type='CD19' then output jnc_CD19_mean_&sysparm.;
    if cell_type='CD4' then output jnc_CD4_mean_&sysparm.;
    if cell_type='CD8' then output jnc_CD8_mean_&sysparm.;
    drop _TYPE_ _FREQ_ cell_type;
run;

data jnc_CD19_mean_2_&sysparm.;
   set jnc_CD19_mean_&sysparm.;
   rename depth=mean_depth_CD19;
run;

data jnc_CD4_mean_2_&sysparm.;
   set jnc_CD4_mean_&sysparm.;
   rename depth=mean_depth_CD4;
run;

data jnc_CD8_mean_2_&sysparm.;
   set jnc_CD8_mean_&sysparm.;
   rename depth=mean_depth_CD8;
run;

data jnc_means_merge_&sysparm.;
   merge jnc_CD19_mean_2_&sysparm. jnc_CD4_mean_2_&sysparm. jnc_CD8_mean_2_&sysparm.;
   by event_id;
run;

/* making permenant */

data mysas.splicing_means_by_celltype_&sysparm.;
set jnc_means_merge_&sysparm.;
run;

