/***** IMPORT EXONS AND FORMAT *****/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import exons */

    data WORK.EXONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exons.csv'
delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chrom $4. ;
       informat start best32. ;
       informat stop best32. ;
       informat strand $1. ;
       informat exon_id $39. ;
       informat transcript_id $1522. ;
       informat gene_id $36. ;
       format chrom $4. ;
       format start best12. ;
       format stop best12. ;
       format strand $1. ;
       format exon_id $39. ;
       format transcript_id $1522. ;
       format gene_id $36. ;
    input
                chrom $
                start
                stop
                strand $
                exon_id $
                transcript_id $
                gene_id $
    ;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  run;
   
/* Make permenant set - can use this later for things */


data av.hg19_aceview_exons;
   set exons;
   exon_length=stop-start;
run;

/* format exons and calculate length */

data exons_formatted;
   set exons;
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

/* test this with PINK1 and some other genes */


data pink1_exons;
  set exons_formatted;
  if gene_id='PINK1' or gene_id='PTPN22' or gene_id='GCH1' or gene_id='UBASH3A' or gene_id='tutee';
run;


data exon_grouping;
   retain exon_group region_start region_stop;
   set pink1_exons;
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
           else do; * exon is inside than current region's coordinates;
                region_start=region_start;
                region_stop=region_stop;
                end;
           exon_group=exon_group;
           end;
       else do; * exon does not overlap;
           exon_group=exon_group+1;
           region_start=start;
           region_stop=stop;
           end;
       end;
run;

/* test works, run on full exon set */

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


/* Seems to work! */

/* Flag longest exon in group as reference exon */
/* If two are the longest, then make the first of the two the reference exon */

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
   keep gene_id exon_group start stop;
   rename start=ref_start stop=ref_stop;
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

/* Add in alternative donor/acceptor flags */

data exon_flag_alt_donor_acceptor;
   set exon_groups_w_refs;
   if start=ref_start then flag_alt_acceptor=0;
   else flag_alt_acceptor=1;
   if stop=ref_stop then flag_alt_donor=0;
   else flag_alt_donor=1;
run;

/* Checks to make sure that no reference exons have an alt donor or acceptor */
proc freq data=exon_flag_alt_donor_acceptor;
   tables flag_reference_exon*flag_alt_donor;
run;


/*

  flag_reference_exon
            flag_alt_donor

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 100277 | 228823 | 329100
           |  14.78 |  33.72 |  48.49
           |  30.47 |  69.53 |
           |  22.29 | 100.00 |
  ---------+--------+--------+
         1 | 349564 |      0 | 349564
           |  51.51 |   0.00 |  51.51
           | 100.00 |   0.00 |
           |  77.71 |   0.00 |
  ---------+--------+--------+
  Total      449841   228823   678664
              66.28    33.72   100.00

*/
proc freq data=exon_flag_alt_donor_acceptor;
   tables flag_reference_exon*flag_alt_acceptor;
run;

/*
 flag_reference_exon
           flag_alt_acceptor

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  99105 | 229995 | 329100
          |  14.60 |  33.89 |  48.49
          |  30.11 |  69.89 |
          |  22.09 | 100.00 |
 ---------+--------+--------+
        1 | 349564 |      0 | 349564
          |  51.51 |   0.00 |  51.51
          | 100.00 |   0.00 |
          |  77.91 |   0.00 |
 ---------+--------+--------+
 Total      448669   229995   678664
             66.11    33.89   100.00

*/

/* looks good! */

/* Make donor exon set */

data donor_exons;
   set exon_flag_alt_donor_acceptor;
   keep exon_group start stop exon_id gene_id flag_short_exon flag_alt_donor;
   rename exon_group=donor_group
          start=donor_exon_start
          stop=donor_exon_stop
          exon_id=donor_exon
          gene_id=donor_gene
          flag_short_exon=flag_short_donor;
run;

/* Make acceptor exon set */

data acceptor_exons;
   set exon_flag_alt_donor_acceptor;
   keep exon_group start stop exon_id gene_id flag_short_exon flag_alt_acceptor;
   rename exon_group=acceptor_group
          start=acceptor_exon_start
          stop=acceptor_exon_stop
          exon_id=acceptor_exon
          gene_id=acceptor_gene
          flag_short_exon=flag_short_acceptor;
run;

/* Make permenant */

data splice.donor_exons;
   set donor_exons;
run;

data splice.acceptor_exons;
  set acceptor_exons;
run;

data splice.exon_list;
    set exon_flag_alt_donor_acceptor;
run;

/* remove temp datasets */

proc datasets noprint;
  delete EXONS exons_formatted pink1_exons exon_grouping
  set_exon_reference  exon_flag_alt_donor_acceptor
  donor_exons acceptor_exons;
run;
quit;


