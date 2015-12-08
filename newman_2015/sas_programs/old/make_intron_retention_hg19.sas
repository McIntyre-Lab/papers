/* Creating intron retention events from exons */

/* Load libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/media/jrbnewman/SAS_WRK1/';

/* We only need the first and last exon from each group for this. Since the reference exon is the first in group, this means we just need to get the last exon per ref_exon group*/

/* sort by ref_exon chr ref_exon_start ref_exon_stop exon_id start stop */

proc sort data=splice2.exon_info_hg19;
   by ref_exon chr ref_exon_start ref_exon_stop exon_id start stop;
run;

data exon_group_info;
   set splice2.exon_info_hg19;
   by ref_exon;
   if last.ref_exon then output;
run;

/* Making dataset a bit smaller (removing unneeded variables) */
/* Ref_exon and exon_id info define the retention event */
/* stop is the latest stop pos in the group, ref_exon_start is the earlier start pos in group*/ 

data retention_exongroup_info;
   retain ref_exon;
   retain exon_id;
   retain gene_cat;
   retain chr;
   retain ref_exon_start;
   retain stop;
   retain strand;
   set exon_group_info;
   rename
       ref_exon=start_exon
       exon_id=stop_exon
       stop=event_stop
       ref_exon_start=event_start;
   drop start ref_exon_stop exon_number ref_exon_strand;
run;


/* sort split file on + or - strand */

proc sort data=retention_exongroup_info;
   by gene_cat chr event_start event_stop;
run;

data retention_exonsgroup_plus  retention_exonsgroup_minus oops;
   set retention_exongroup_info;
   if strand="+" then output retention_exonsgroup_plus;
   else if strand="-" then output retention_exonsgroup_minus;
   else output oops; *0 obs!;
run;

/* set to flags for last exon_group per gene */

data exonsgroup_plus_flag_last;
   set retention_exonsgroup_plus;
   by gene_cat;
   if last.gene_cat then flag_last_group=1;
   else flag_last_group=0;
run;

data exonsgroup_minus_flag_last;
   set retention_exonsgroup_minus;
   by gene_cat;
   if first.gene_cat then flag_last_group=1;
   else flag_last_group=0;
run;


/* Calculate intron size and add start and stop positions per event */
/* if fusion or intron size <38bp then take full length of feature */

/* Need to sort descending for plus */
/* Minus is fine as is already in reverse order */

/* For antisense fusions */
data exonsgroup_minus_length;
   set exonsgroup_minus_flag_last;

   intron_start=lag(event_stop)+1; *adding 1 otherwise length will include 1bp of previous fusion;
   intron_size=event_start-intron_start; *get intron length;
   fusion_size=event_stop-event_start; *get fusion length;

   if intron_size gt 37 then intron_ret_stop = event_start-37; *set event start pos;
   else intron_ret_stop=event_start-intron_size;
   if fusion_size gt 38 then intron_ret_start = event_start+38; *set event stop pos;
   else intron_ret_start=event_start+fusion_size;
   if intron_size lt 0 then flag_skip_intron=1;
   else if flag_last_group=1 then flag_skip_intron=1;
   else flag_skip_intron=0;
run;


/* For sense fusions */
proc sort data=exonsgroup_plus_flag_last;
   by gene_cat chr descending event_start;
run;

data exonsgroup_plus_length;
   set exonsgroup_plus_flag_last;
   intron_stop=lag(event_start)-1; *subtracting 1 otherwise length will include 1bp of next fusion;
   intron_size=intron_stop-event_start; *get intron length;
   fusion_size=event_stop-event_start; *get fusion length;

   if intron_size gt 38 then intron_ret_stop = event_stop+38; *set event start pos;
   else intron_ret_stop=event_stop+intron_size;
   if fusion_size gt 37 then intron_ret_start = event_stop-37; *set event stop pos;
   else intron_ret_start=event_stop-fusion_size;

   if intron_size lt 0 then flag_skip_intron=1;
   else if flag_last_group=1 then flag_skip_intron=1;
   else flag_skip_intron=0;
run;

/* check that intron_ret events are not bigger than 150bp */

/* data size_check;
   set exonsgroup_minus_length;
   ret_event_length=intron_ret_stop-intron_ret_start;
   if flag_skip_intron=0;
run; */

/* Drop all but necessary variables */
/* Drop exon groups with no retention event */

data exongroups_plus_cat;
   set exonsgroup_plus_length;
   if flag_skip_intron=0;
   drop flag_skip_intron fusion_size intron_size intron_stop flag_last_group;
run;

data exongroups_minus_cat;
   set exonsgroup_minus_length;
   if flag_skip_intron=0;
   drop flag_skip_intron fusion_size intron_size intron_start flag_last_group;
run;


/* Concatenate datasets */
/* Add intron retention id */

data intron_retention_all;
   length event_id $100.;
   set exongroups_plus_cat exongroups_minus_cat;
   if strand="+" then event_id=catx('|',stop_exon,'intron');
   else event_id=catx('|',start_exon,'intron');
   flag_intron_retention=1;
run;


/* Make data permenant */

data splice2.intron_retention_events_hg19;
    set intron_retention_all;
run;

/* data intron_retention_events2;
   set intron_retention_all;
   if strand="+" then do;
      bed_start=intron_ret_start;
      bed_stop=intron_ret_stop;
      end;
   else do;
      bed_start=intron_ret_stop;
      bed_stop=intron_ret_start;
      end;
run;


data intron_retention_bed;
   retain chr;
   retain bed_start;
   retain bed_stop;
   retain event_id;
   length event_id2 $100.;
   set intron_retention_events2;
   event_id2=tranwrd(event_id,"|","_");
   keep chr bed_start bed_stop event_id2;
run;

proc export data=intron_retention_bed
	outfile='/home/jrbnewman/McLab/junction_annotations/generated_files/dmel_r557_catalogue_intron_retention.bed'
	dbms=tab replace;
        putnames=no;
	run;

*/





