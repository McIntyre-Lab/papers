/* Code to check if in eachb exon group there is only one reference exon */

libname splice '/mnt/data/splice/';
libname av '/home/jrbnewman/ugfi_share/SHARE/McIntyre_Lab/useful_human_data/aceview_hg19/sas_data';

/* format exons and calculate length */

data exons_formatted;
   set av.hg19_aceview_exons;
   exon_length=stop-start;
   if exon_length lt 37 then flag_short_exon=1;
   else flag_short_exon=0;
   drop transcript_id;
run;


/* Group exons */
/* Sort exons by start and end
   Calculate length
   Group exons (see IR script for logic)
   Flag longest
   If longest, reference donor/acceptor
   If not longest, alternative donor/acceptor */


proc sort data=exons_formatted;
   by gene_id chrom start stop;
run;



data exon_grouping;
   retain exon_group region_start region_stop;
   set exons_formatted;
   by gene_id;
   if first.gene_id then do;
       exon_group=1;
       region_start=start;
       region_stop=stop;
       end;
   else do;
       if start le region_stop then do; *check to see if exon overlaps with previous exon;
           * if exon overlaps previous see which has the bigger coordinates;
           if stop ge region_stop then do; *expand exon region;
                region_start=region_start;
                region_stop=stop;
                end;
           else do; * exon is inside the current region's coordinates;
                region_start=region_start;
                region_stop=region_stop;
                end;
           exon_group=exon_group;
           end;
       else do; * exon does not overlap with exon region in memory, then begin a new exon region;
           exon_group=exon_group+1;
           region_start=start;
           region_stop=stop;
           end;
       end;
run;

/* Sort exon groups by descending exon length and start position */

proc sort data=exon_grouping;
   by gene_id exon_group descending exon_length start;
run;


data set_exon_reference;
   set exon_grouping;
   by gene_id exon_group;
   if first.exon_group then flag_reference_exon=1;
   else flag_reference_exon=0;
run;

/* Need to get references for each group */

proc sort data=set_exon_reference;
   by gene_id exon_group region_start region_stop;
run;

data exon_references;
   set set_exon_reference;
   by gene_id exon_group;
   if flag_reference_exon=1;
   keep gene_id exon_group start stop exon_length;
   rename start=ref_start stop=ref_stop exon_length=ref_length;
run;

proc sort data=exon_references;
   by gene_id exon_group;
run;


proc sort data=set_exon_reference;
  by gene_id exon_group;
run;

data exon_groups_w_refs oops1 oops2;
   merge exon_references (in=in1) set_exon_reference (in=in2);
   by gene_id exon_group;
   if in1 and in2 then output exon_groups_w_refs;
   else if in1 then output oops1;
   else output oops2;
run;

/* Flag if exon length is the same as the reference length */


data flag_exon_length;
  set exon_groups_w_refs;
  if exon_length=ref_length then flag_length=1;
  else flag_length=0;
run;


/* Keep exons that are the same length as their reference exon */


data flag_exon_length2;
 set flag_exon_length;
 if flag_length=1;
run;

/* Get frequency of each exon group. These should all be 1! */

proc sort data=flag_exon_length2;
   by gene_id exon_group;
run;

proc freq data=flag_exon_length2 noprint;
   by gene_id;
   tables exon_group / out=exon_length_check;
run;

proc sort data=exon_length_check;
  by descending count;
run;

/* 37 non-ref exons are the same length as the reference exon */

