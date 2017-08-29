/*Need to sum technical replicates across samples for analysis. 

Example of data:
NAME         SUBJECT_ID      CELL_TYPE
PC_0001         M001            CD19
PC_0002         M001            CD8
PC_0003         M001            CD4

Sum by Name and transcript_id to get each transcript with each subject and cell type
*/

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';



*Test on one transcript;

*data single;
*  set mysas.all_counts_w_key;
*  if transcript_id = "ENST00000207157";
*  run;


*Sort by transcript_id;
proc sort data=mysas.all_counts_w_key;
  by Name fusion_id;
  run;

*Sum technical reps for each sample;
proc means data = mysas.all_counts_w_key noprint;
  *class Name;
  by Name fusion_id;
  var region_depth;
  id region_length;
  output out=mysas.mapped_reads_summed_byevent sum=depth;
run;

* Recalculate APN and make dataset permanent;
data mysas.counts_by_event;
   set mysas.mapped_reads_summed_byevent;
   if name ne " " ; 
   apn = depth / region_length;
   run;

* calculate total mapped reads by subject;

data mysas.mapped_reads ;
  set mysas.all_counts_w_key;
  keep subject_id mapped_reads; 
  run;

proc sort data=mysas.mapped_reads nodup;
  by subject_id mapped_reads;
  run;

proc means data=mysas.mapped_reads;
  by subject_id;
  var mapped_reads;
  output out=mysas.total_mapped_reads_sum sum=total_mapped_reads;
  run;


*export for check;
proc export data=mysas.mapped_reads_summed_byevent
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/mapped_reads_summed.csv'
    label dbms=csv replace;
    run;

proc export data=mysas.total_mapped_reads_sum
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/total_mapped_reads.csv'
    label dbms=csv replace;
    run;


