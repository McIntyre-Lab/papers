libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* I want to look at only gene pairs where each gene is ONLY in that pair.,

E.g. if gene X is in the pairs X|Y and X|Z, I don't want to use gene X */

data gene_list;
  set event.multigene_frags_on;
  keep gene_id;
run;

proc sort data=gene_list nodup;
  by gene_id;
run;


data gene_list_stack;
  length unstack_gene_id $15.;
  set gene_list;
  do i=1 by 1 while(scan(gene_id,i,"|") ^= '');
     unstack_gene_id=scan(gene_id,i,"|");
     output;
     end;
run;

proc freq data=gene_list_stack noprint;
  tables unstack_gene_id / out=count_pairs;
run;

data genes_to_drop;
  set count_pairs;
  where count > 1;
  keep unstack_gene_id;
run;

proc sort data=genes_to_drop;
  by unstack_gene_id;
proc sort data=gene_list_stack;
  by unstack_gene_id;
run;

data event.pairs_to_keep;
  merge gene_list_stack (in=in1) genes_to_drop (in=in2);
  by unstack_gene_id;
  if in2 then delete;
  keep gene_id;
run;

