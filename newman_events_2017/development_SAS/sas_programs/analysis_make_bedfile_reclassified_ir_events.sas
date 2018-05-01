ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Create and export a BED file of reclassified IR events */

data ir_events;
  set event.ir_reclassification_v2;
  where flag_low_expressed=0;
  if flag_possible_unprocessed=1 
  or flag_possible_novel_donor=1
  or flag_possible_ir=1 then output;
  keep event_id;
run;


data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data ir2gene;
   set evspl.splicing_events_annot_refseq;
   where flag_intron_retention=1;
   keep event_id gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop event_size
        feature1_size feature2_size;
run;

proc sort data=ir2gene;
   by gene_id;
proc sort data=exp_genes;
   by gene_id;
run;

data ir2gene_exp;
  merge exp_genes (in=in1) ir2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=ir2gene_exp;
  by event_id;
proc sort data=ir_events;
  by event_id;
run;

data ir2gene_exp_on;
  merge ir2gene_exp (in=in1) ir_events (in=in2);
  by event_id;
  if in1 and in2;
run; *34780 unannotated events detected;



/* Assemble and export BED file */


data info_for_bed;
   length score $1.;
   length color $7.;
   set ir2gene_exp_on;
   score='.';
   color='255,0,0';
   block1_start=0;
      num_blocks=1;
      if strand='+' then do;
         block1_length=event_size+1;
         totalstart1=feature1_start;
         totalstop1=feature2_stop+1;
         totalstart2=feature1_start;
         totalstop2=feature2_stop+1;
         end;
      if strand='-' then do;
         block1_length=event_size+1;
         totalstart1=feature1_start-1;
         totalstop1=feature2_stop;
         totalstart2=feature1_start-1;
         totalstop2=feature2_stop;
         end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand
        score color num_blocks block1_length  block1_start ;
run;


/* BED12 format 


chrom           chr
totalStart      feature1_start
totalStop       feature2_stop
name            event_id
score           "."
strand          strand
totalStart      feature1_start
totalStop       feature2_stop
color           255,0,0
num             "1" if intron retention, "2" if junction
lengths         block lengths (total if intron_retention, feature1_stop-feature1_start,feature2_stop-feature2_start
starts          if retention then 0, if junction then (0,feature2_start-feature1_start)

*/


data assemble_bed;
    retain chr;
    retain totalstart1;
    retain totalstop1;
    retain event_id;
    retain score;
    retain strand;
    retain totalstart2;
    retain totalstop2;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
        lengths_cat=put(block1_length, 3.);
        starts_cat=put(block1_start, 1.);
    drop block1_length block1_start;
run;


proc sort data=assemble_bed;
   by chr totalstart1 totalstop1 strand;
run;



/* output as BED */

proc export data=assemble_bed
        outfile='!MCLAB/event_analysis/references/reclassified_ir_events_npc.bed'
        dbms=tab replace;
        putnames=no;
        run;

