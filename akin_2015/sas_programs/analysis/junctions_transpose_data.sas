/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* Transpose depth by sample */

data jnc_data_for_tpose;
   set mysas.jnc_counts_w_flags_&sysparm.;
   keep subject fusion_id region_depth;
run;

proc sort data=jnc_data_for_tpose;
    by fusion_id subject;
run;

proc transpose data=jnc_data_for_tpose out=mysas.jnc_data_flip_&sysparm.;
   by fusion_id;
   var region_depth;
   id subject;
run;



