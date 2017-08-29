/* Merge counts with flags and drop "off" events */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data2';

data splicing_flags;
   set mysas.counts_by_splicing_w_flags_&sysparm.;
   keep event_id flag_CD19_on flag_CD4_on flag_CD8_on;
run;

data design_file;
    set mysas.design_file_se;
    keep Name cell_type subject_id;
run;

proc sort data=splicing_flags nodups;
   by event_id;
proc sort data=mysas.splicing_counts_w_zeros_&sysparm.;
    by event_id;
run;

/* Merge flags with counts, but only keep counts with flags! */

data splicing_counts_w_flags;
   merge splicing_flags (in=in1) mysas.splicing_counts_w_zeros_&sysparm.;
   by event_id;
   if in1;
run;

proc sort data=design_file;
   by Name;
run;

proc sort data=splicing_counts_w_flags;
   by Name;
run;


/* Merge flags with counts, but only keep counts with data! */

data splicing_counts_flags_info;
   merge design_file (in=in1) splicing_counts_w_flags (in=in2);
   by Name;
   if in1 and in2 then output;
run;

/* Remove events not "on" in cell type where observation is also cell type */


data splicing_counts_pruned;
    set splicing_counts_flags_info;
    if cell_type='CD4' and flag_cd4_on=0 then delete;
    if cell_type='CD8' and flag_cd8_on=0 then delete;
    if cell_type='CD19' and flag_cd19_on=0 then delete;
run;


/* Make permenant */


data mysas.splicing_counts_for_eqtls_&sysparm.;
   set splicing_counts_pruned;
run;

