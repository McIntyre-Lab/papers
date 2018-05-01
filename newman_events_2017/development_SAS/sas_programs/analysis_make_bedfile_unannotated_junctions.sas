ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Create and export a BED file of detected unannotated junctions */

data unannot_on;
   set event.splicing_on_apn_gt0;
   where flag_splicing_on=1;
   keep event_id;
run;

data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data unannot2gene;
   set evspl.splicing_events_annot_refseq;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   keep event_id gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop event_size
        feature1_size feature2_size;
run;

proc sort data=unannot2gene;
   by gene_id;
proc sort data=exp_genes;
   by gene_id;
run;

data unannot2gene_exp;
  merge exp_genes (in=in1) unannot2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=unannot2gene_exp;
  by event_id;
proc sort data=unannot_on;
  by event_id;
run;

data unannot2gene_exp_on;
  merge unannot2gene_exp (in=in1) unannot_on (in=in2);
  by event_id;
  if in1 and in2;
run; *15796 unannotated events detected;

/* Make permenant */

data event.unannot_events_dtctd_exp_genes;
   set unannot2gene_exp_on;
run;

/* Assemble and export BED file */

data info_for_bed;
   length score $1.;
   length color $7.;
   set unannot2gene_exp_on;
   score='.';
   color='255,0,0';
   block1_start=0;
      num_blocks=2;
          totalstart1=feature1_start;
          totalstop1=feature2_stop;
          totalstart2=feature1_start;
          totalstop2=feature2_stop;
      block1_length=feature1_size;
      block2_length=feature2_size;
      block2_start=feature2_start-feature1_start;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
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
    lengths_cat=catx(',', block1_length, block2_length);
    starts_cat=catx(',', block1_start, block2_start);
    drop block1_length block2_length block1_start block2_start;
run;

proc sort data=assemble_bed;
   by chr totalstart1 totalstop1 strand;
run;


/* output as BED */

proc export data=assemble_bed
        outfile='!MCLAB/event_analysis/references/unannotated_junctions_detected_npc.bed'
        dbms=tab replace;
        putnames=no;
        run;

