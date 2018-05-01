ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Subset IKZF3, FAM129C, FLT4 and IL8 (three diff, one not diff) */

data subset_data;
  set eventloc.t1d_exons_mee_rank_cell_subj_v2;
  where gene_id in ("IKZF3","FAM129C","IL8","FLT4");
run;

/* Make permenant */

data event.t1d_genes_for_mee_testing;
   set subset_data;
run;

   
