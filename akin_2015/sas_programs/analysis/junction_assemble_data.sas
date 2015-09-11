/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* Assemble data! */

/* get annotations */

data jnc_annot;
   set mysas.splicing_events_annotations;
   keep event_id event_type gene_id flag_junction_annotated flag_intron_retention num_skipped_exons flag_exonskip flag_alt_donor flag_alt_acceptor;
   rename event_id=fusion_id;
run;


data junc_info;
   set mysas.jnc_counts_w_flags_&sysparm.;
   keep fusion_id region_length flag_control_on flag_treat_on flag_all_on;
run;

proc sort data=junc_info nodup;
   by fusion_id;
run;

/* pieces:  jnc_annot    
/* build summarized data (pre-ANOVA) */

proc sort data=jnc_annot;
    by fusion_id;
run;

proc sort data=mysas.jnc_data_flip_&sysparm.;
    by fusion_id;
run;

proc sort data=mysas.jnc_means_merge_&sysparm.;
    by fusion_id;
run;

proc sort data=mysas.junc_counts_merge_&sysparm.;
    by fusion_id;
run;


data junc_counts_depth oops;
   merge mysas.junc_counts_merge_&sysparm. (in=in1) mysas.jnc_data_flip_&sysparm. (in=in2);
   by fusion_id;
   if in1 and in2 then output junc_counts_depth;
   else if in1 then output oops;
   else do;
       num_samples_exp=0;
       num_samples_present=0;
       output junc_counts_depth;
      end;
run;

data junc_counts_flags oops1 oops2;
   merge junc_info (in=in1) junc_counts_depth (in=in2);
   by fusion_id;
   if in1 and in2 then output junc_counts_flags;
    else if in1 then output oops1;
    else output oops2;
run;

data junc_counts_means oops1 oops2;
   merge mysas.jnc_means_merge_&sysparm. (in=in1) junc_counts_flags (in=in2);
   by fusion_id;
   if in1 and in2 then output junc_counts_means;
    else if in1 then output oops1;
    else output oops2;
run;

data junc_counts_w_annot oops1 oops2;
   merge jnc_annot (in=in1) junc_counts_means (in=in2);
   by fusion_id;
   if in1 and in2 then output junc_counts_w_annot;
    else if in1 then output oops1;
    else output oops2;
run;


/*flag goi make perm */

data mysas.junc_counts_w_annot_&sysparm.;
   set junc_counts_w_annot;
      if gene_id='PAX6'
      or gene_id='FGFR2'
      or gene_id='CD44'
      or gene_id='CTNND1'
      or gene_id='TMX2andC11orf31andCTNND1'
      or gene_id='ENAH'
      or gene_id='SLC37A2'
      or gene_id='ARHGEF11'
      or gene_id='FOXJ3'
      or gene_id='FAM50A'
      or gene_id='PSENEN'
      or gene_id='PSENENandLIN37'
      or gene_id='ECT2'
      or gene_id='NCSTN'
      or gene_id='SLC1A2'
      or gene_id='NCRNA00085'
      or gene_id='NCRNA00077'
      or gene_id='HAS2AS' then flag_goi=1;
      else flag_goi=0;
run;


