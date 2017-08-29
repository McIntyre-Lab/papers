/*** IMPORT JUNCTION DATA ***/

/* import libraries */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';


/********************* FLAG COUNTS *************************/


/*data splicing_counts;
   set mysas.all_counts_w_key;
   length sample_id $33.;
   index1=scan(sample_id,1,'_');
   index2=scan(sample_id,2,'_');
   index3=scan(sample_id,4,'_');
   index4=scan(sample_id,5,'_');
   index5=scan(sample_id,6,'_');
   sample_id2=catx('_',index1,index2,index3,index4,index5);
   drop index1 index2 index3 index4 index5 sample_id library index barcode Name pool repeat_pool sample_description cell_type subject_id lane flowcell total_reads;
   rename sample_id2=sample_id;
run;
*/


proc sort data=mysas.splicing_counts;
   by sample_id;
run;

proc sort data=mysas.design_file;
   by sample_id;
run;

data mysas.splicing_counts_w_key oops1 oops2;
   merge mysas.splicing_counts (in=in1) mysas.design_file (in=in2);
   by sample_id;
   if in1 and in2 then output mysas.splicing_counts_w_key;
   else if in1 then output oops1;
   else output oops2;
run;


/* Sum tech reps */

proc sort data=mysas.splicing_counts_w_key;
   by Event_id Name;
run;

proc means data=mysas.splicing_counts_w_key noprint;
   by Name event_id;
   var region_depth;
   output out=mysas.splicing_counts_sum sum=;
run;

